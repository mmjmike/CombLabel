import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", "-f", help='Specify file for which to calculate stats', type=str)
    args = parser.parse_args()
    try:
        f = open(args.file, "r")
    except FileNotFoundError:
        print("File '{}' doesn't exist".format(args.file))
        return
    cycle read lines until eof
        if






if __name__ == "__main__":
    main()