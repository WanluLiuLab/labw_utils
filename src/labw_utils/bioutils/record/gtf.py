from typing import List, Optional

from labw_utils.bioutils.record.feature import Feature, DEFAULT_GTF_QUOTE_OPTIONS, VALID_GTF_QUOTE_OPTIONS, \
    feature_repr
from labw_utils.commonutils.str_utils import to_dict


def format_string(
        feature: Feature,
        quote: str = DEFAULT_GTF_QUOTE_OPTIONS
):
    if quote not in VALID_GTF_QUOTE_OPTIONS:
        raise ValueError(f"Invalid quoting option {quote}, should be one in {VALID_GTF_QUOTE_OPTIONS}.")
    attribute_full_str = ""
    for k, v in zip(feature.attribute_keys, feature.attribute_values):
        attr_str = feature_repr(v)
        if quote == "blank":
            if any(map(lambda blank_or_sep: blank_or_sep in attr_str, " \t\n\f\r;")):
                attr_str = f"\"{attr_str}\""
        elif quote == "string" and isinstance(v, str):
            attr_str = f"\"{attr_str}\""
        elif quote == "all":
            attr_str = f"\"{attr_str}\""
        attribute_full_str = f"{attribute_full_str}{k} " + attr_str + "; "
    return ("\t".join((
        feature_repr(feature.seqname),
        feature_repr(feature.source),
        feature_repr(feature.feature),
        feature_repr(feature.start),
        feature_repr(feature.end),
        feature_repr(feature.score),
        "." if feature.strand is None else ("+" if feature.strand else "-"),
        feature_repr(feature.frame),
        attribute_full_str[:-1]
    )))


def parse_record(
        in_str: str,
        skip_fields: Optional[List[str]] = None,
        included_attributes: Optional[List[str]] = None,
) -> Feature:
    if skip_fields is None:
        skip_fields = []
    line_split = in_str.rstrip('\n\r').split('\t')
    required_fields = line_split[0:-1]
    attributes = to_dict(
        line_split[-1],
        field_sep=' ',
        record_sep=';',
        quotation_mark='\"\'',
        resolve_str=True
    )

    if included_attributes is not None:
        attributes = {k: v for k, v in attributes.items() if k in included_attributes}

    return Feature(
        seqname=required_fields[0],
        source=required_fields[1] if required_fields[1] != "." and "source" not in skip_fields else None,
        feature=required_fields[2] if required_fields[2] != "." and "feature" not in skip_fields else None,
        start=int(required_fields[3]),
        end=int(required_fields[4]),
        score=float(required_fields[5]) if required_fields[5] != "." and "score" not in skip_fields else None,
        strand=required_fields[6] == "+" if required_fields[6] != "." and "strand" not in skip_fields else None,
        frame=int(required_fields[7]) if required_fields[7] != "." and "frame" not in skip_fields else None,
        attribute=attributes
    )
