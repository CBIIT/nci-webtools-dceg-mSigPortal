import os
import rpy2


def testPython(args):
    print(args.params)
    return {True, args.params}


def testR(args):


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    subparsers = parser.add_subparsers()

    testPythonParser = subparsers.add_parser('testPython')
    testPythonParser.add_argument('params')
    testPythonParser.set_defaults(func=testPython)

    args = parser.parse_args()
    args.func(args)

# parser.add_argument('-f', '--function', help='specify function to call')
# args = parser.parse_args()
# if (args.debug):
#     config = Util('config.dev.ini')
# else:
