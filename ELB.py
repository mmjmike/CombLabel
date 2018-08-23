#!/usr/bin/python3

import argparse







def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", "-c", help='Specify the config file', type=str)
    args = parser.parse_args()



if __name__ == "__main__":
    main()
