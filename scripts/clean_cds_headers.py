import argparse
import re

def main(input, output):
    with open(input, 'r') as s_in, open(output, 'w') as s_out:
        for line in s_in:
            if '>' in line:
                protein_id = re.search("(?<=cds_).{5,}?(?=_)", line).group(0)
                s_out.write(''.join(['>',protein_id,'\n']))
            else:
                s_out.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    args = parser.parse_args()
    main(**vars(args))
