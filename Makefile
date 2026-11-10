.PHONY: install py-install web-install py-lint py-typecheck py-test py-ci web-lint web-typecheck web-test web-build web-ci ci fmt dev-web contract contract-check

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

# ---- Everything -----------------------------------------------------------

ci: py-ci contract-check web-ci

fmt:
	uv run ruff format .
	uv run ruff check --fix .

dev-web:
	cd web && npm run dev
