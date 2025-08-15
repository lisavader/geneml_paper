import argparse
import re

def main(input, output):
    pseudo_count = 0
    with open(input, 'r') as s_in, open(output, 'w') as s_out:
        for line in s_in:
            if '>' in line:
                if "pseudo=true" in line:
                    pseudo_count += 1
                    continue
                try:
                    protein_id = re.search("(?<=cds_)[A-Z0-9_.]{6,}?(?=_)", line).group(0)
                except AttributeError:
                    raise ValueError(f"Invalid line: {line}")
                s_out.write(''.join(['>',protein_id,'\n']))
            else:
                s_out.write(line)
    print(f"Done.\nSkipped {pseudo_count} pseudogenes.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    args = parser.parse_args()
    main(**vars(args))
