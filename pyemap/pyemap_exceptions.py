__all__ = [
    'PyeMapException', 'PyeMapUserResidueException', 'PyeMapGraphException', 'PyeMapShortestPathException',
    'PyeMapMiningException', 'PyeMapParseException', 'PyeMapGraphDatabaseException'
]


class PyeMapException(Exception):
    pass


class PyeMapUserResidueException(PyeMapException):
    pass


class PyeMapGraphException(PyeMapException):
    pass


class PyeMapShortestPathException(PyeMapException):
    pass


class PyeMapMiningException(PyeMapException):
    pass


class PyeMapParseException(PyeMapException):
    pass


class PyeMapGraphDatabaseException(PyeMapException):
    pass
