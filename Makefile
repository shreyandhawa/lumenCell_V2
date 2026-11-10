.PHONY: install py-install web-install py-lint py-typecheck py-test py-ci web-lint web-typecheck web-test web-build web-ci ci fmt dev-web

install: py-install web-install

py-install:
	uv sync

web-install:
	cd web && npm ci

# ---- Python ---------------------------------------------------------------

py-lint:
	uv run ruff check .
	uv run ruff format --check .

py-typecheck:
	uv run mypy

py-test:
	uv run pytest -m "not real_data and not slow"

py-ci: py-lint py-typecheck py-test

# ---- Web ------------------------------------------------------------------

web-lint:
	cd web && npm run -s lint

web-typecheck:
	cd web && npm run -s typecheck

web-test:
	cd web && npm run -s test

web-build:
	cd web && npm run -s build

web-ci: web-lint web-typecheck web-test web-build

# ---- Everything -----------------------------------------------------------

ci: py-ci web-ci

fmt:
	uv run ruff format .
	uv run ruff check --fix .

dev-web:
	cd web && npm run dev
