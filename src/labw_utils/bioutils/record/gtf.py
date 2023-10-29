"""
``labw_utils.bioutils.record.gtf`` -- GTF parsers

.. versionadded:: 1.0.2
"""

__all__ = ("format_string", "parse_record")

from labw_utils.bioutils.record.feature import (
    FeatureInterface,
    DEFAULT_GTF_QUOTE_OPTIONS,
    VALID_GTF_QUOTE_OPTIONS,
    feature_repr,
    GtfAttributeValueType,
    strand_repr,
    Feature,
)
from labw_utils.commonutils.str_utils import to_dict
from labw_utils.typing_importer import List, Optional


class GTFParsingError(ValueError):
    """
    General GTF parsing errors.

    .. versionadded:: 1.0.2"""

    pass


def _format_attribute_str(attribute_key: str, attribute_value: GtfAttributeValueType, quote: str) -> str:
    if isinstance(attribute_value, List):
        return "; ".join(
            _format_attribute_str(attribute_key, _single_value, quote) for _single_value in attribute_value
        )
    attr_str = feature_repr(attribute_value).replace(";", "!").replace('"', "!").replace("'", "!")
    if quote == "blank":
        if any(map(lambda blank_or_sep: blank_or_sep in attr_str, " \t\n\f\r;")):
            attr_str = f'"{attr_str}"'
    elif quote == "string" and isinstance(attribute_value, str):
        attr_str = f'"{attr_str}"'
    elif quote == "all":
        attr_str = f'"{attr_str}"'
    return f"{attribute_key} {attr_str}"


def format_string(feature: FeatureInterface, quote: str = DEFAULT_GTF_QUOTE_OPTIONS):
    """
    Format :py:class:`Feature` to GTF.

    :param feature: Feature to be formatted.
    :param quote: Quoting policy.

    :raise ValueError: On invalid quoting options.

    .. versionadded:: 1.0.2
    """
    if quote not in VALID_GTF_QUOTE_OPTIONS:
        raise ValueError(f"Invalid quoting option {quote}, should be one in {VALID_GTF_QUOTE_OPTIONS}.")
    attribute_full_str = (
        "; ".join(
            _format_attribute_str(attribute_key, attribute_value, quote)
            for attribute_key, attribute_value in zip(feature.attribute_keys, feature.attribute_values)
        )
        + ";"
    )
    return "\t".join(
        (
            feature_repr(feature.seqname),
            feature_repr(feature.source),
            feature_repr(feature.feature),
            feature_repr(feature.start),
            feature_repr(feature.end),
            feature_repr(feature.score),
            strand_repr(feature.strand),
            feature_repr(feature.frame),
            attribute_full_str,
        )
    )


def parse_record(
    in_str: str,
    skip_fields: Optional[List[str]] = None,
    included_attributes: Optional[List[str]] = None,
) -> FeatureInterface:
    """
    Parse record string to :py:class:`Feature`.

    :param in_str: String to be parsed.
    :param skip_fields: Explicitly skip optional features to reduce space.
    :param included_attributes: Explicitly include attributes to reduce space. Other attributes are discarded.

    :raises GTFParsingError: On invalid record.

    .. versionadded:: 1.0.2
    """
    if skip_fields is None:
        skip_fields = []
    line_split = in_str.rstrip("\n\r").split("\t")
    if len(line_split) != 9:
        raise GTFParsingError(f"Illegal GTF record '{in_str}': Should have 9 fields, here only {len(line_split)}")
    required_fields = line_split[0:-1]
    attributes = to_dict(line_split[-1], field_sep=" ", record_sep=";", quotation_mark="\"'", resolve_str=True)

    if included_attributes is not None:
        attributes = {k: v for k, v in attributes.items() if k in included_attributes}

    try:
        return Feature(
            seqname=required_fields[0],
            source=required_fields[1] if required_fields[1] != "." and "source" not in skip_fields else None,
            feature=required_fields[2] if required_fields[2] != "." and "feature" not in skip_fields else None,
            start=int(required_fields[3]),
            end=int(required_fields[4]),
            score=float(required_fields[5]) if required_fields[5] != "." and "score" not in skip_fields else None,
            strand=required_fields[6] == "+" if required_fields[6] != "." and "strand" not in skip_fields else None,
            frame=int(required_fields[7]) if required_fields[7] != "." and "frame" not in skip_fields else None,
            attribute=attributes,
        )
    except ValueError as e:
        raise GTFParsingError(f"Illegal GTF record '{in_str}'") from e
