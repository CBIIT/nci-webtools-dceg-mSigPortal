import os
import rpy2
import json
from python import test


def pythonHandler(args):
    req = json.loads(args.params)
    fn = getattr(test, req['fn'])
    args = req['args']
    print(json.dumps({'return': fn(args)}), end='')


def rHandler(args):
    print('rHander: ', args.params)


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    subparsers = parser.add_subparsers()

    pythonParser = subparsers.add_parser('pythonHandler')
    pythonParser.add_argument('params')
    pythonParser.set_defaults(func=pythonHandler)

    rParser = subparsers.add_parser('rHandler')
    rParser.add_argument('params')
    rParser.set_defaults(func=rHandler)

    args = parser.parse_args()
    args.func(args)
