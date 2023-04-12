"""
labw_utils.stdlib_helper.argparse_helper -- Argument parser with enhanced help formatter

Following is an example using "normal" formatter:

>>> parser = argparse.ArgumentParser(prog="prog", description="description")
>>> _ = parser.add_argument("p", type=int, help="p-value")
>>> _ = parser.add_argument("-o", required=True, type=str, help="output filename", default="/dev/stdout")
>>> _ = parser.add_argument("--flag", action="store_true", help="flag")
>>> print(parser.format_help())
usage: prog [-h] -o O [--flag] p
<BLANKLINE>
description
<BLANKLINE>
positional arguments:
  p           p-value
<BLANKLINE>
optional arguments:
  -h, --help  show this help message and exit
  -o O        output filename
  --flag      flag
<BLANKLINE>


Following is an example using enhanced formatter:

>>> parser = ArgumentParserWithEnhancedFormatHelp(prog="prog", description="description")
>>> _ = parser.add_argument("p", type=int, help="p-value")
>>> _ = parser.add_argument("-o", required=True, type=str, help="output filename", default="/dev/stdout")
>>> _ = parser.add_argument("--flag", action="store_true", help="flag")
>>> print(parser.format_help())
description
<BLANKLINE>
SYNOPSIS: prog [-h] -o O [--flag] p
<BLANKLINE>
PARAMETERS:
  p
              [REQUIRED] Type: int;
              p-value
<BLANKLINE>
OPTIONS:
  -h, --help
              [OPTIONAL]
              show this help message and exit
  -o O
              [REQUIRED] Type: str; Default: /dev/stdout
              output filename
  --flag
              [OPTIONAL] Default: False
              flag
<BLANKLINE>
"""

__all__ = (
    "ArgumentParserWithEnhancedFormatHelp",
)

import argparse


class _EnhancedHelpFormatter(argparse.HelpFormatter):
    def _expand_help(self, action: argparse.Action):
        params = {**vars(action), "prog": self._prog}
        for name in list(params):
            if params[name] is argparse.SUPPRESS:
                del params[name]
        for name in list(params):
            if hasattr(params[name], '__name__'):
                params[name] = params[name].__name__
        if params.get('choices') is not None:
            choices_str = ', '.join([str(c) for c in params['choices']])
            params['choices'] = choices_str

        help_str = action.help

        if action.required:
            req_opt_prefix = "[REQUIRED] "
        else:
            req_opt_prefix = "[OPTIONAL] "

        if action.type is None:
            dtype_prefix = ""
        else:
            dtype_prefix = "Type: " + action.type.__name__ + "; "

        default_prefix = ""
        if '%(default)' not in action.help:
            if action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    default_prefix = f'Default: {action.default} '
        return (req_opt_prefix + dtype_prefix + default_prefix).strip() + "\n" + help_str % params

    def _format_action(self, action):
        help_position = min(
            self._action_max_length + 2,
            self._max_help_position
        )
        action_header = '%*s' % (self._current_indent, '') + self._format_action_invocation(action) + "\n"
        parts = [action_header]
        if action.help:
            help_text = self._expand_help(action)
            help_lines = help_text.split("\n")
            parts.append('%*s%s\n' % (help_position, '', help_lines[0]))
            for line in help_lines[1:]:
                parts.append('%*s%s\n' % (help_position, '', line))
        elif not action_header.endswith('\n'):
            parts.append('\n')
        for subaction in self._iter_indented_subactions(action):
            parts.append(self._format_action(subaction))
        return self._join_parts(parts)


_ACTION_GROUP_TILE_REPLACEMENT_DICT = {
    "optional arguments": "OPTIONS",  # Python < 3.10
    "options": "OPTIONS",  # Python >= 3.10
    "positional arguments": "PARAMETERS"
}


class ArgumentParserWithEnhancedFormatHelp(argparse.ArgumentParser):
    def format_help(self) -> str:
        formatter = _EnhancedHelpFormatter(prog=self.prog)
        formatter.add_text(self.description)

        formatter.add_usage(
            usage=self.usage,
            actions=self._actions,
            groups=self._mutually_exclusive_groups,
            prefix="SYNOPSIS: "
        )

        for action_group in self._action_groups:
            formatter.start_section(_ACTION_GROUP_TILE_REPLACEMENT_DICT.get(action_group.title, action_group.title))
            formatter.add_arguments(action_group._group_actions)
            formatter.end_section()

        formatter.add_text(self.epilog)
        return formatter.format_help()
