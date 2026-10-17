#!/usr/bin/env python3
"""Validate MGF metadata and create a DeepMASS2-compatible copy (Python 3.8+)."""
import argparse
import hashlib
import json
import math
import re
import sys
from pathlib import Path

VERSION = "1.1.0"
UNSAFE = re.compile(r'[<>:"/\\|?*\x00-\x1f]')
RESERVED = re.compile(r'^(CON|PRN|AUX|NUL|COM[1-9]|LPT[1-9])(?:\.|$)', re.I)


def positive_mass(value, label):
    try:
        mass = float(value)
    except (TypeError, ValueError):
        raise ValueError("Invalid {}: {!r}".format(label, value))
    if not math.isfinite(mass) or mass <= 0:
        raise ValueError("{} must be finite and positive".format(label))
    return mass


def parse_mgf(text):
    """Return blocks with original lines, metadata, and exact peak lines."""
    blocks, current = [], None
    for number, line in enumerate(text.splitlines(), 1):
        marker = line.strip().upper()
        if marker == 'BEGIN IONS':
            if current is not None:
                raise ValueError("Nested BEGIN IONS at line {}".format(number))
            current = {'start': number - 1, 'fields': {}, 'peaks': [], 'field_lines': {}}
        elif marker == 'END IONS':
            if current is None:
                raise ValueError("Unexpected END IONS at line {}".format(number))
            if not current['peaks']:
                raise ValueError("Empty spectrum ending at line {}".format(number))
            current['end'] = number - 1
            blocks.append(current)
            current = None
        elif current is not None:
            if not line.strip() or line.lstrip().startswith(('#', ';', '!')):
                continue
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip().upper()
                if key in current['fields']:
                    raise ValueError("Duplicate {} at line {}".format(key, number))
                current['fields'][key] = value.strip()
                current['field_lines'][key] = number - 1
            else:
                columns = line.split()
                if len(columns) < 2:
                    raise ValueError("Invalid peak at line {}".format(number))
                positive_mass(columns[0], 'peak m/z at line {}'.format(number))
                try:
                    intensity = float(columns[1])
                except ValueError:
                    raise ValueError("Invalid intensity at line {}".format(number))
                if not math.isfinite(intensity) or intensity < 0:
                    raise ValueError("Invalid intensity at line {}".format(number))
                current['peaks'].append(line)
    if current is not None:
        raise ValueError('Missing END IONS')
    if not blocks:
        raise ValueError('No spectra found')
    return blocks


def safe_name(name):
    if (not name or UNSAFE.search(name) or name.endswith(('.', ' '))
            or RESERVED.match(name) or len(name) > 180):
        raise ValueError('Unsafe or excessively long output name: {!r}'.format(name))
    return name


def peak_digest(blocks):
    # Preserve spectrum boundaries and original peak text in this fingerprint.
    payload = json.dumps([b['peaks'] for b in blocks], ensure_ascii=False).encode('utf-8')
    return hashlib.sha256(payload).hexdigest()


def convert(raw, default_ionmode=None):
    text = raw.decode('utf-8-sig')
    newline = '\r\n' if b'\r\n' in raw else '\n'
    lines = text.splitlines()
    blocks = parse_mgf(text)
    names, modes, additions, replacements = set(), {}, {}, {}
    changes = {'precursor_mz_added': 0, 'compound_name_added': 0, 'compound_name_replaced': 0,
               'ionmode_added': 0, 'ionmode_normalized': 0}
    for index, block in enumerate(blocks, 1):
        fields = block['fields']
        try:
            extra = []
            pepmass = fields.get('PEPMASS', '').split()
            explicit = fields.get('PRECURSOR_MZ', '')
            if not pepmass and not explicit:
                raise ValueError('Missing PEPMASS and PRECURSOR_MZ')
            mass_text = pepmass[0] if pepmass else explicit
            mass = positive_mass(mass_text, 'precursor m/z')
            if explicit and abs(positive_mass(explicit, 'PRECURSOR_MZ') - mass) > 1e-8:
                raise ValueError('Conflicting PEPMASS and PRECURSOR_MZ')
            if not explicit:
                if 'PRECURSOR_MZ' in fields:
                    raise ValueError('Empty PRECURSOR_MZ')
                extra.append('PRECURSOR_MZ=' + mass_text)
                changes['precursor_mz_added'] += 1
            ionmode = fields.get('IONMODE', '').lower()
            if not ionmode:
                if not default_ionmode:
                    raise ValueError('Missing IONMODE; supply --ionmode explicitly')
                ionmode = default_ionmode
                if 'IONMODE' in fields:
                    replacements[block['field_lines']['IONMODE']] = 'IONMODE=' + ionmode
                else:
                    extra.append('IONMODE=' + ionmode)
                changes['ionmode_added'] += 1
            elif ionmode not in ('positive', 'negative'):
                raise ValueError('IONMODE must be positive or negative')
            else:
                replacements[block['field_lines']['IONMODE']] = 'IONMODE=' + ionmode
                if fields['IONMODE'] != ionmode:
                    changes['ionmode_normalized'] += 1
            modes[ionmode] = modes.get(ionmode, 0) + 1
            # Use the original feature identifier verbatim, without a prefix or renumbering.
            name = fields.get('FEATURE_ID')
            if not name:
                raise ValueError('Missing FEATURE_ID; cannot assign COMPOUND_NAME')
            if 'COMPOUND_NAME' not in fields:
                extra.append('COMPOUND_NAME=' + name)
                changes['compound_name_added'] += 1
            elif fields['COMPOUND_NAME'] != name:
                replacements[block['field_lines']['COMPOUND_NAME']] = 'COMPOUND_NAME=' + name
                changes['compound_name_replaced'] += 1
            safe_name(name)
            if name.casefold() in names:
                raise ValueError('Duplicate output name (case-insensitive): ' + name)
            names.add(name.casefold())
            additions[block['end']] = extra
        except ValueError as error:
            raise ValueError('Spectrum {}: {}'.format(index, error)) from error
    output = []
    for index, line in enumerate(lines):
        output.extend(additions.get(index, []))
        output.append(replacements.get(index, line))
    encoded = (newline.join(output) + newline).encode('utf-8')
    new_blocks = parse_mgf(encoded.decode('utf-8'))
    fingerprint = peak_digest(blocks)
    if fingerprint != peak_digest(new_blocks):
        raise ValueError('Internal error: peak data changed')
    report = {'tool_version': VERSION, 'naming_rule': 'COMPOUND_NAME=FEATURE_ID', 'spectra': len(blocks),
              'peak_lines': sum(len(b['peaks']) for b in blocks),
              'unique_names': len(names), 'ionmode_counts': modes, 'changes': changes,
              'input_sha256': hashlib.sha256(raw).hexdigest(),
              'output_sha256': hashlib.sha256(encoded).hexdigest(),
              'peak_data_sha256': fingerprint, 'all_peak_lines_unchanged': True}
    return encoded, report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path, help='UTF-8 MGF input')
    parser.add_argument('output', type=Path, nargs='?', help='New output file; never overwritten')
    parser.add_argument('--ionmode', choices=['positive', 'negative'],
                        help='Fill missing IONMODE only; never override existing polarity')
    parser.add_argument('--check', action='store_true', help='Validate and preview changes without writing')
    parser.add_argument('--report', type=Path, help='Write JSON audit report to a new file')
    parser.add_argument('--version', action='version', version=VERSION)
    args = parser.parse_args(argv)
    if not args.check and args.output is None:
        parser.error('output is required unless --check is used')
    if args.check and args.report:
        parser.error('--check does not write files; omit --report')
    try:
        if not args.check:
            for target in (args.output, args.report):
                if target is not None and target.exists():
                    raise ValueError('Refusing to overwrite: ' + str(target))
            if args.report is not None and args.report.absolute() == args.output.absolute():
                raise ValueError('Output and report must be different files')
        encoded, report = convert(args.source.read_bytes(), args.ionmode)
        # Basenames keep the audit report suitable for public sharing.
        report.update({'input_file': args.source.name,
                       'output_file': args.output.name if args.output else None,
                       'check_only': args.check})
        if not args.check:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            with args.output.open('xb') as handle:
                handle.write(encoded)
            if args.report:
                args.report.parent.mkdir(parents=True, exist_ok=True)
                with args.report.open('x', encoding='utf-8') as handle:
                    json.dump(report, handle, ensure_ascii=False, indent=2)
                    handle.write('\n')
        print(json.dumps(report, ensure_ascii=False, indent=2))
        return 0
    except (OSError, ValueError) as error:
        print('Error: {}'.format(error), file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
