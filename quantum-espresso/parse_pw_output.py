#!/usr/bin/env python


from os import path


# -- the information we want to extract from the output file of PROFESS
KPTS_PATTERN = ')   bands (ev):'
KPTS_KEYWORD = 'kpoints'
ENERGY_PATTERN = '!    total energy              ='
ENERGY_KEYWORD = 'energy'
TIME_PATTERN = '     PWSCF        :'
TIME_KEYWORD = 'time'

pattern_to_keyword = {
    KPTS_PATTERN: KPTS_KEYWORD,
    ENERGY_PATTERN: ENERGY_KEYWORD,
    TIME_PATTERN: TIME_KEYWORD,
}

pattern_l = [
    TIME_PATTERN,
    KPTS_PATTERN,
    ENERGY_PATTERN,
]

extract_fct_l = [
    lambda line: line.split('CPU')[0].split(':')[-1].strip(),
    lambda line: [float(item) for item in line.split()],
    lambda line: float(line.split()[-2]),
]


# -- Core extraction routine
def get_data_from_patterns(
    lines,
    pattern_l,
    extraction_fct_l,
    ):

    assert len(pattern_l) == len(extraction_fct_l)

    data_dict = {}
    for pattern in pattern_l:
        data_dict[pattern] = []
    for idx, line in enumerate(lines):
        for pattern, extract_fct in zip(pattern_l, extraction_fct_l):
            # -- special case of k points for which the proper line is the
            # second next line underneath
            if pattern in line and pattern == KPTS_PATTERN:
                try:
                    data = extract_fct(lines[idx + 2])
                    data_dict[pattern] += [data]
                except:
                    print "could not extract data"
                    print "from line : '%s'" % (line,)
                    continue
            elif pattern in line and pattern != KPTS_PATTERN:
                try:
                    data = extract_fct(line)
                    data_dict[pattern] += [data]
                except:
                    print "could not extract data"
                    print "from line : '%s'" % (line,)
                    continue
    return data_dict


# -- Front-end output parsing routine
def parse_pw_output(
    output_fname,
    ):

    assert path.isfile(output_fname)

    with open(output_fname, 'r') as fin:
        lines = fin.readlines()

    data_dict = get_data_from_patterns(
            lines,
            pattern_l,
            extract_fct_l,
            )
    final_dict = {}
    for pattern, data in data_dict.items():
        final_dict[pattern_to_keyword[pattern]] = data

    return final_dict
