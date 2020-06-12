def test(args):
    print(args.param)
    return {True, args.param}


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    subparsers = parser.add_subparsers()

    testParser = subparsers.add_parser('test')
    testParser.add_argument('param')
    testParser.set_defaults(func=test)

    args = parser.parse_args()
    args.func(args)

# parser.add_argument('-f', '--function', help='specify function to call')
# args = parser.parse_args()
# if (args.debug):
#     config = Util('config.dev.ini')
# else:
