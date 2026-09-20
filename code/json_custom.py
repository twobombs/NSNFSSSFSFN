import json
import abc


class json_custom_serializable(abc.ABC):
    """
    Inherit from this class, and define the member functions as_dict()
    and to_dict(), in order to get json serialization.
    """

    # the sage.rings.integer.Integer thing is only a convenient special
    # case so that sage doesn't pester me with long ints not being
    # serializable to json
    encoders = {
            'sage.rings.integer.Integer': int,
            'sage.rings.finite_rings.integer_mod.IntegerMod_int': int,
            'sage.rings.finite_rings.integer_mod.IntegerMod_gmp': int,
            'sage.rings.real_double_element_gsl.RealDoubleElement_gsl': float
            }
    decoders = {}

    def __init_subclass__(cls):
        key = f"{cls.__module__}.{cls.__name__}"
        cls.encoders[key] = cls.as_dict
        cls.decoders[key] = cls.from_dict

    @abc.abstractmethod
    def as_dict(self):
        """
        return a dictionary representation of self, so that
        type(self).to_dict(self.as_dict()) is equivalents to self (which
        is weaker that requiring that all members be equal -- as it turns
        out, in most cases we don't need the full memory.
        """
        pass

    @classmethod
    @abc.abstractmethod
    def from_dict(cls):
        """
        convert the dictionary representation {obj} into an object of the
        current class
        """
        pass


# idea from https://mathspp.com/blog/custom-json-encoder-and-decoder

class MyEncoder(json.JSONEncoder):
    def default(self, obj):
        module = type(obj).__module__
        name = type(obj).__name__
        key = f"{module}.{name}"
        if m := json_custom_serializable.encoders.get(key):
            d = m(obj)
            # The returned serialization should really be a dictionary,
            # except in the case where we would get something that can be
            # converted both ways, like sage.rings.integer.Integer vs int
            if type(d) is dict:
                d["__extended_json_type__"] = key
            return d
        else:
            super().default(obj)


class MyDecoder(json.JSONDecoder):
    def __init__(self, **kwargs):
        kwargs["object_hook"] = self.object_hook
        super().__init__(**kwargs)

    def object_hook(self, obj):
        try:
            name = obj["__extended_json_type__"]
            return json_custom_serializable.decoders[name](obj)
        except (KeyError, AttributeError):
            return obj


def dumps(*args, **kwargs):
    return json.dumps(*args, cls=MyEncoder, **kwargs)


def loads(*args, **kwargs):
    return json.loads(*args, cls=MyDecoder, **kwargs)


def dump(*args, **kwargs):
    return json.dump(*args, cls=MyEncoder, **kwargs)


def load(*args, **kwargs):
    return json.load(*args, cls=MyDecoder, **kwargs)
