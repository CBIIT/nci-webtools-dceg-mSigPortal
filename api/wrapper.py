import os
import json
from rpy2.robjects import r

from python import test


def pythonHandler(args):
    req = json.loads(args.params)
    fn = getattr(test, req['fn'])
    args = req['args']
    result = fn(args)
    print(json.dumps({'return': result}), end='')


def rHandler(args):
    req = json.loads(args.params)
    fn = req['fn']
    args = req['args']
    r.source('api/R/test.R')
    result = str(r[fn](args))
    print(json.dumps({'return': result}), end='')


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
