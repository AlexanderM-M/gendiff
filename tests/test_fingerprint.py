from semantiseq.fingerprint import digest


def test_digest_is_typed_and_canonical() -> None:
    assert digest({"b": 2, "a": [1, None]}) == digest({"a": (1, None), "b": 2})
    assert digest([1]) != digest(["1"])
    assert digest(True) != digest(1)
