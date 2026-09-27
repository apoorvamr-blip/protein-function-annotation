"""HTTP client retry/backoff/rate-limit tests with a fake requests session."""

from __future__ import annotations

import pytest
import requests

from protein_annotator.utils.http import HttpClient, HttpError, RateLimiter, ServiceUnavailableError


class FakeSession:
    def __init__(self, outcomes):
        self.outcomes = list(outcomes)
        self.headers: dict[str, str] = {}
        self.requests = []

    def request(self, method, url, **kwargs):
        self.requests.append((method, url, kwargs))
        outcome = self.outcomes.pop(0)
        if isinstance(outcome, Exception):
            raise outcome
        return outcome

    def close(self):
        pass


def response(status, text="", headers=None):
    r = requests.Response()
    r.status_code = status
    r._content = text.encode()
    r.headers.update(headers or {})
    return r


def client(session, sleeps, max_retries=3):
    return HttpClient(user_agent="test-agent/1.0", max_retries=max_retries, backoff_factor=1.0,
                      session=session, sleep=sleeps.append)


def test_sets_user_agent():
    session = FakeSession([response(200, "ok")])
    client(session, []).get("https://x.example/a")
    assert session.headers["User-Agent"] == "test-agent/1.0"


def test_retries_transient_status_then_succeeds():
    sleeps: list[float] = []
    session = FakeSession([response(503), response(429, headers={"Retry-After": "7"}), response(200, "ok")])
    assert client(session, sleeps).get("https://x.example/a").text == "ok"
    assert len(session.requests) == 3
    assert sleeps[1] == 7.0  # Retry-After honoured
    assert 1.0 <= sleeps[0] <= 1.25  # backoff_factor * 2**0 plus jitter


def test_retries_connection_errors_with_exponential_backoff():
    sleeps: list[float] = []
    session = FakeSession([requests.ConnectionError("boom")] * 3 + [response(200, "ok")])
    client(session, sleeps).get("https://x.example/a")
    assert sleeps[0] < sleeps[1] < sleeps[2]


def test_gives_up_after_max_retries():
    session = FakeSession([response(502)] * 3)
    with pytest.raises(ServiceUnavailableError):
        client(session, [], max_retries=2).get("https://x.example/a")
    assert len(session.requests) == 3


def test_client_errors_are_not_retried():
    session = FakeSession([response(404, "missing")])
    with pytest.raises(HttpError) as info:
        client(session, []).get("https://x.example/a")
    assert info.value.status_code == 404
    assert len(session.requests) == 1


def test_rate_limiter_spaces_requests_per_host():
    now = [0.0]
    sleeps: list[float] = []

    def sleep(seconds):
        sleeps.append(seconds)
        now[0] += seconds

    limiter = RateLimiter({"a.example": 10.0}, clock=lambda: now[0], sleep=sleep)
    limiter.wait("a.example")
    now[0] += 3
    limiter.wait("a.example")
    limiter.wait("b.example")  # no interval configured
    assert sleeps == [pytest.approx(7.0)]
