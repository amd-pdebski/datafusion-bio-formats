//! Predicate pushdown helpers for [`crate::table_provider::VcfTableProvider`].
//! Only a fixed set of core columns is advertised as pushable; evaluation must match that set
//! so comet-bio / DataFusion can safely use `df.filter` without incorrect results.

use datafusion::arrow::datatypes::Schema;
use datafusion::common::ScalarValue;
use datafusion::logical_expr::{Between, Expr, Operator, expr::InList};
use noodles::vcf::{Header, Record};
use noodles::vcf::variant::record::{AlternateBases, Filters, Ids};
use noodles::vcf::variant::Record as VariantRecord;
use std::sync::Arc;

/// `end` column value consistent with [`crate::physical_exec`] batch building.
pub(crate) fn vcf_end_column_value(record: &Record, header: &Header) -> u32 {
    vcf_end_column_value_dyn(record as &dyn VariantRecord, header)
}

pub(crate) fn vcf_end_column_value_dyn(record: &dyn VariantRecord, header: &Header) -> u32 {
    let ref_len = record.reference_bases().len();
    let alt_len = record.alternate_bases().len();
    if ref_len == 1
        && alt_len == 1
        && record
            .reference_bases()
            .iter()
            .map(|c| c.unwrap())
            .all(|c| c == b'A' || c == b'C' || c == b'G' || c == b'T')
        && record
            .alternate_bases()
            .iter()
            .map(|c| c.unwrap())
            .all(|c| c.eq("A") || c.eq("C") || c.eq("G") || c.eq("T"))
    {
        record.variant_start().unwrap().unwrap().get() as u32
    } else {
        record.variant_end(header).unwrap().get() as u32
    }
}

/// Core VCF columns we both accept for pushdown and evaluate during scan (INFO / FORMAT excluded).
const PUSHABLE_COLUMNS: &[&str] =
    &["chrom", "start", "end", "qual", "filter", "id", "ref", "alt"];

fn is_pushable_column(name: &str) -> bool {
    PUSHABLE_COLUMNS.iter().any(|&c| c == name)
}

pub fn can_push_down_filter(expr: &Expr, schema: &Arc<Schema>) -> bool {
    match expr {
        Expr::BinaryExpr(binary_expr) if matches!(binary_expr.op, Operator::And) => {
            can_push_down_filter(&binary_expr.left, schema)
                && can_push_down_filter(&binary_expr.right, schema)
        }
        Expr::BinaryExpr(binary_expr) => can_push_down_binary_expr(binary_expr, schema),
        Expr::Between(between_expr) => can_push_down_between_expr(between_expr, schema),
        Expr::InList(in_list_expr) => can_push_down_in_list_expr(in_list_expr, schema),
        _ => false,
    }
}

/// Evaluates filter expressions against a parsed VCF record.
pub fn evaluate_filters_against_record(
    record: &Record,
    header: &Header,
    filters: &[Expr],
) -> bool {
    if filters.is_empty() {
        return true;
    }
    filters
        .iter()
        .all(|filter| evaluate_single_filter(record, header, filter))
}

fn evaluate_single_filter(record: &Record, header: &Header, filter: &Expr) -> bool {
    match filter {
        Expr::BinaryExpr(binary_expr) if matches!(binary_expr.op, Operator::And) => {
            evaluate_single_filter(record, header, &binary_expr.left)
                && evaluate_single_filter(record, header, &binary_expr.right)
        }
        Expr::BinaryExpr(binary_expr) => evaluate_binary_filter(record, header, binary_expr),
        Expr::Between(between_expr) => evaluate_between_filter(record, header, between_expr),
        Expr::InList(in_list_expr) => evaluate_in_list_filter(record, header, in_list_expr),
        _ => true,
    }
}

fn record_end_u64(record: &Record, header: &Header) -> u64 {
    u64::from(vcf_end_column_value(record, header))
}

fn evaluate_binary_filter(
    record: &Record,
    header: &Header,
    binary_expr: &datafusion::logical_expr::BinaryExpr,
) -> bool {
    if let Expr::Column(column) = &*binary_expr.left {
        if let Expr::Literal(literal, _) = &*binary_expr.right {
            let field_name = column.name.as_str();
            return match field_name {
                "chrom" => evaluate_string_comparison(
                    record.reference_sequence_name(),
                    literal,
                    &binary_expr.op,
                ),
                "start" => match record.variant_start() {
                    Some(Ok(pos)) => {
                        let v = pos.get() as f64;
                        evaluate_numeric_comparison(v, literal, &binary_expr.op)
                    }
                    _ => false,
                },
                "end" => {
                    let v = record_end_u64(record, header) as f64;
                    evaluate_numeric_comparison(v, literal, &binary_expr.op)
                }
                "qual" => match record.quality_score() {
                    Some(Ok(q)) => {
                        evaluate_numeric_comparison(f64::from(q), literal, &binary_expr.op)
                    }
                    _ => evaluate_numeric_comparison(0.0, literal, &binary_expr.op),
                },
                "filter" => {
                    let filters_str = record
                        .filters()
                        .iter(header)
                        .map(|v| v.unwrap_or(".").to_string())
                        .collect::<Vec<String>>()
                        .join(";");
                    evaluate_string_comparison(&filters_str, literal, &binary_expr.op)
                }
                "id" => {
                    let id = record
                        .ids()
                        .iter()
                        .map(|v| v.to_string())
                        .collect::<Vec<String>>()
                        .join(";");
                    evaluate_string_comparison(&id, literal, &binary_expr.op)
                }
                "ref" => evaluate_string_comparison(record.reference_bases(), literal, &binary_expr.op),
                "alt" => {
                    let alt = record
                        .alternate_bases()
                        .iter()
                        .map(|v| v.unwrap_or(".").to_string())
                        .collect::<Vec<String>>()
                        .join("|");
                    evaluate_string_comparison(&alt, literal, &binary_expr.op)
                }
                _ => false,
            };
        }
    }
    false
}

fn evaluate_between_filter(
    record: &Record,
    header: &Header,
    between_expr: &Between,
) -> bool {
    if let Expr::Column(column) = &*between_expr.expr {
        if let (Expr::Literal(low_literal, _), Expr::Literal(high_literal, _)) =
            (&*between_expr.low, &*between_expr.high)
        {
            let field_name = column.name.as_str();
            let negated = between_expr.negated;
            return match field_name {
                "start" => match record.variant_start() {
                    Some(Ok(pos)) => evaluate_between_comparison(
                        pos.get() as f64,
                        low_literal,
                        high_literal,
                        negated,
                    ),
                    _ => false,
                },
                "end" => {
                    let v = record_end_u64(record, header) as f64;
                    evaluate_between_comparison(v, low_literal, high_literal, negated)
                }
                "qual" => match record.quality_score() {
                    Some(Ok(q)) => evaluate_between_comparison(
                        f64::from(q),
                        low_literal,
                        high_literal,
                        negated,
                    ),
                    _ => false,
                },
                _ => false,
            };
        }
    }
    false
}

fn evaluate_in_list_filter(record: &Record, header: &Header, in_list_expr: &InList) -> bool {
    if let Expr::Column(column) = &*in_list_expr.expr {
        let field_name = column.name.as_str();
        let values: Vec<String> = in_list_expr
            .list
            .iter()
            .filter_map(|e| match e {
                Expr::Literal(sv, _) => match sv {
                    ScalarValue::Utf8(Some(s)) => Some(s.clone()),
                    ScalarValue::LargeUtf8(Some(s)) => Some(s.clone()),
                    ScalarValue::UInt32(Some(v)) => Some(v.to_string()),
                    ScalarValue::Int32(Some(v)) => Some(v.to_string()),
                    ScalarValue::Int64(Some(v)) => Some(v.to_string()),
                    ScalarValue::Float32(Some(v)) => Some(v.to_string()),
                    ScalarValue::Float64(Some(v)) => Some(v.to_string()),
                    _ => None,
                },
                _ => None,
            })
            .collect();

        let contains = match field_name {
            "chrom" => values.contains(&record.reference_sequence_name().to_string()),
            "filter" => {
                let filters_str = record
                    .filters()
                    .iter(header)
                    .map(|v| v.unwrap_or(".").to_string())
                    .collect::<Vec<String>>()
                    .join(";");
                values.contains(&filters_str)
            }
            "id" => {
                let id = record
                    .ids()
                    .iter()
                    .map(|v| v.to_string())
                    .collect::<Vec<String>>()
                    .join(";");
                values.contains(&id)
            }
            "ref" => values.contains(&record.reference_bases().to_string()),
            "alt" => {
                let alt = record
                    .alternate_bases()
                    .iter()
                    .map(|v| v.unwrap_or(".").to_string())
                    .collect::<Vec<String>>()
                    .join("|");
                values.contains(&alt)
            }
            _ => false,
        };
        return if in_list_expr.negated {
            !contains
        } else {
            contains
        };
    }
    false
}

fn evaluate_string_comparison(
    record_value: &str,
    literal: &ScalarValue,
    op: &Operator,
) -> bool {
    let literal_value = match literal {
        ScalarValue::Utf8(Some(s)) => s.as_str(),
        ScalarValue::LargeUtf8(Some(s)) => s.as_str(),
        _ => return false,
    };
    match op {
        Operator::Eq => record_value == literal_value,
        Operator::NotEq => record_value != literal_value,
        _ => false,
    }
}

fn evaluate_numeric_comparison(record_value: f64, literal: &ScalarValue, op: &Operator) -> bool {
    let literal_value = match literal {
        ScalarValue::UInt32(Some(val)) => *val as f64,
        ScalarValue::Float32(Some(val)) => f64::from(*val),
        ScalarValue::Float64(Some(val)) => *val,
        ScalarValue::Int32(Some(val)) => *val as f64,
        ScalarValue::Int64(Some(val)) => *val as f64,
        _ => return false,
    };
    match op {
        Operator::Eq => (record_value - literal_value).abs() < f64::EPSILON,
        Operator::NotEq => (record_value - literal_value).abs() >= f64::EPSILON,
        Operator::Lt => record_value < literal_value,
        Operator::LtEq => record_value <= literal_value,
        Operator::Gt => record_value > literal_value,
        Operator::GtEq => record_value >= literal_value,
        _ => false,
    }
}

fn evaluate_between_comparison(
    record_value: f64,
    low_literal: &ScalarValue,
    high_literal: &ScalarValue,
    negated: bool,
) -> bool {
    let low_value = match low_literal {
        ScalarValue::UInt32(Some(val)) => *val as f64,
        ScalarValue::Float32(Some(val)) => f64::from(*val),
        ScalarValue::Float64(Some(val)) => *val,
        ScalarValue::Int32(Some(val)) => *val as f64,
        ScalarValue::Int64(Some(val)) => *val as f64,
        _ => return false,
    };
    let high_value = match high_literal {
        ScalarValue::UInt32(Some(val)) => *val as f64,
        ScalarValue::Float32(Some(val)) => f64::from(*val),
        ScalarValue::Float64(Some(val)) => *val,
        ScalarValue::Int32(Some(val)) => *val as f64,
        ScalarValue::Int64(Some(val)) => *val as f64,
        _ => return false,
    };
    let between = record_value >= low_value && record_value <= high_value;
    if negated {
        !between
    } else {
        between
    }
}

fn can_push_down_binary_expr(
    binary_expr: &datafusion::logical_expr::BinaryExpr,
    schema: &Arc<Schema>,
) -> bool {
    if let Expr::Column(column) = &*binary_expr.left {
        let field_name = column.name.as_str();
        if !is_pushable_column(field_name) {
            return false;
        }
        if schema.field_with_name(field_name).is_err() {
            return false;
        }
        if !matches!(&*binary_expr.right, Expr::Literal(..)) {
            return false;
        }
        return match field_name {
            "chrom" | "filter" | "id" | "ref" | "alt" => {
                matches!(binary_expr.op, Operator::Eq | Operator::NotEq)
            }
            "start" | "end" | "qual" => matches!(
                binary_expr.op,
                Operator::Eq
                    | Operator::NotEq
                    | Operator::Lt
                    | Operator::LtEq
                    | Operator::Gt
                    | Operator::GtEq
            ),
            _ => false,
        };
    }
    false
}

fn can_push_down_between_expr(between_expr: &Between, schema: &Arc<Schema>) -> bool {
    if let Expr::Column(column) = &*between_expr.expr {
        let field_name = column.name.as_str();
        if !is_pushable_column(field_name) {
            return false;
        }
        if !matches!(field_name, "start" | "end" | "qual") {
            return false;
        }
        if schema.field_with_name(field_name).is_err() {
            return false;
        }
        return matches!(&*between_expr.low, Expr::Literal(..))
            && matches!(&*between_expr.high, Expr::Literal(..));
    }
    false
}

fn can_push_down_in_list_expr(in_list_expr: &InList, schema: &Arc<Schema>) -> bool {
    if let Expr::Column(column) = &*in_list_expr.expr {
        let field_name = column.name.as_str();
        if !is_pushable_column(field_name) {
            return false;
        }
        if !matches!(field_name, "chrom" | "filter" | "id" | "ref" | "alt") {
            return false;
        }
        if schema.field_with_name(field_name).is_err() {
            return false;
        }
        return in_list_expr
            .list
            .iter()
            .all(|expr| matches!(expr, Expr::Literal(..)));
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::datatypes::{DataType, Field};
    use datafusion::logical_expr::{col, lit};

    #[test]
    fn pushdown_whitelist_rejects_info_column() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("AC", DataType::Int32, true),
        ]));
        let expr = col("AC").eq(lit(1_i32));
        assert!(!can_push_down_filter(&expr, &schema));
    }

    #[test]
    fn pushdown_accepts_chrom_eq() {
        let schema = Arc::new(Schema::new(vec![Field::new(
            "chrom",
            DataType::Utf8,
            false,
        )]));
        let expr = col("chrom").eq(lit("chr1"));
        assert!(can_push_down_filter(&expr, &schema));
    }
}
