.PHONY: install py-install py-lint py-typecheck py-test py-ci ci fmt

install: py-install

py-install:
	uv sync

py-lint:
	uv run ruff check .
	uv run ruff format --check .

py-typecheck:
	uv run mypy

py-test:
	uv run pytest -m "not real_data and not slow"

py-ci: py-lint py-typecheck py-test

ci: py-ci

fmt:
	uv run ruff format .
	uv run ruff check --fix .
