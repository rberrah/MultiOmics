import assert from 'node:assert/strict';
import { createRequire } from 'node:module';

const requireFromSveltePlugin = createRequire(import.meta.resolve('@sveltejs/vite-plugin-svelte'));
const merge = requireFromSveltePlugin('deepmerge');

const nested = merge(
  { session: null },
  JSON.parse('{"session":{"__proto__":{"admin":true}}}')
);
assert.equal(Object.getPrototypeOf(nested.session), Object.prototype);
assert.equal(nested.session.admin, undefined);

const nullPrototypeTarget = Object.create(null);
nullPrototypeTarget.debug = false;
const nullPrototypeResult = merge(
  nullPrototypeTarget,
  JSON.parse('{"__proto__":{"admin":true}}')
);
assert.equal(Object.getPrototypeOf(nullPrototypeResult), Object.prototype);
assert.equal(nullPrototypeResult.admin, undefined);
assert.equal(nullPrototypeResult.debug, false);

const ownProtoResult = merge(
  JSON.parse('{"__proto__":{"admin":true},"debug":false}'),
  { verbose: true }
);
assert.equal(ownProtoResult.admin, undefined);
assert.equal(Object.hasOwn(ownProtoResult, '__proto__'), false);
assert.equal(ownProtoResult.verbose, true);

const arrayResult = merge(
  { list: [] },
  JSON.parse('{"list":[{"__proto__":{"admin":true}}]}')
);
assert.equal(Object.getPrototypeOf(arrayResult.list[0]), Object.prototype);
assert.equal(arrayResult.list[0].admin, undefined);
assert.equal({}.admin, undefined);

console.log('deepmerge prototype-pollution regression check passed.');
