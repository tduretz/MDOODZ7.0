.PHONY: test client server


# Development

start:
	docker-compose up -d && docker-compose logs -f -t migration

stop:
	docker-compose stop

ci:
	npm ci --no-audit

check-code-style:
	npm run check-code-style

autoformat:
	npm run autoformat

lint:
	npm run lint

lint-fix:
	npm run lint-fix

audit:
	npm audit --audit-level=critical

build:
	npm run build

dev:
	npm run dev:all

client:
	npm run dev:client

server:
	npm run dev:server

push:
	git push origin HEAD


# Testing

test:
	npm run test $(filter-out $@, $(MAKECMDGOALS))
%: # this allows matching received argument as target name to avoid absent target error


# Deployment

setup: push # TODO fix fail with quotes
	ansible-playbook -i ops/prod.ini ops/setup.yml --extra-vars "git_branch=\"$(shell git branch --show-current)\" release_name=\"$(shell git log -1 --pretty=tformat:'%s\n\nCommit: %h\nAuthor: %an')\"" --ask-vault-pass

deploy: push
	ansible-playbook -i ops/prod.ini ops/deploy.yml --extra-vars "git_branch=\"$(shell git branch --show-current)\" release_name=\"$(shell git log -1 --pretty=tformat:'%s\n\nCommit: %h\nAuthor: %an')\"" --ask-vault-pass
