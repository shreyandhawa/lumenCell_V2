.PHONY: install py-install web-install py-lint py-typecheck py-test py-ci web-lint web-typecheck web-test web-build web-ci web-e2e ci fmt dev-web web-dist dev-server contract contract-check py-test-all bench-synth bench-baseline bench-real

install: py-install web-install

py-install:
	uv sync --all-extras

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

# Everything that runs offline, including slow calibration and benchmark tests.
py-test-all:
	uv run pytest -m "not real_data"

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

# End to end: the real server serving the built app, driven by Playwright (Chromium).
web-e2e: web-build
	cd web && npx playwright test

# ---- API contract ---------------------------------------------------------
# contracts/openapi.json is generated from api/models.py, and the TypeScript
# types are generated from it. `contract` regenerates both; `contract-check`
# fails if either committed file is stale.

contract:
	uv run python -m lumencell.api.app contracts/openapi.json
	cd web && npm run -s gen:api

contract-check:
	@tmp=$$(mktemp -d); \
	uv run python -m lumencell.api.app $$tmp/openapi.json > /dev/null && \
	diff -u contracts/openapi.json $$tmp/openapi.json > /dev/null || { echo "contracts/openapi.json is stale: run make contract"; exit 1; }; \
	(cd web && npx openapi-ts -i ../contracts/openapi.json -o $$tmp/generated -p @hey-api/typescript --no-log-file -s) && \
	diff -ru web/src/api/generated $$tmp/generated > /dev/null || { echo "web/src/api/generated is stale: run make contract"; exit 1; }; \
	rm -rf $$tmp; echo "API contract is up to date"

# ---- Benchmarks -----------------------------------------------------------

BENCH_METHODS ?= paga_dpt,palantir,scfates_ppt

# Full synthetic matrix: reports gates, fails on regressions against the baseline.
bench-synth:
	uv run lumencell bench synth --methods $(BENCH_METHODS)

# Re-measure and replace benchmarks/baseline.json (commit the result deliberately).
bench-baseline:
	uv run lumencell bench synth --methods $(BENCH_METHODS) --write-baseline

# Real datasets (download first: uv run lumencell bench fetch all). A release gate, not CI.
bench-real:
	uv run pytest -m real_data -q
	uv run lumencell bench real

# ---- Everything -----------------------------------------------------------

ci: py-ci contract-check web-ci web-e2e

fmt:
	uv run ruff format .
	uv run ruff check --fix .

dev-web:
	cd web && npm run dev

# Build the web app and bundle it into the Python package, so `lumencell launch` serves it.
web-dist: web-build
	rm -rf src/lumencell/web_dist
	cp -R web/dist src/lumencell/web_dist

# Run the server on a synthetic dataset for development (web: make dev-web in another shell).
dev-server:
	uv run lumencell synth .lumencell/dev.h5ad --store-root
	uv run lumencell launch .lumencell/dev.h5ad --no-browser
