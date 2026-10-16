'use strict';
const { test } = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const { createHash } = require('node:crypto');
const oi = require('../scripts/oi_wiki.cjs');
const headers = require('../scripts/publish_headers.cjs');
const data = Buffer.from('%PDF-1.7\ntest\n');
const digest = `sha256:${createHash('sha256').update(data).digest('hex')}`;
const context = { repo: { owner: 'owner', repo: 'repo' }, sha: 'built',
  ref: 'refs/heads/master', eventName: 'push', payload: { repository: { default_branch: 'master' } } };
const core = { info() {}, warning() {}, setOutput() {} };
const asset = { id: 7, name: 'OI-wiki.pdf', state: 'uploaded', size: data.length, digest };
function fixture(t, bytes = data) {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'xcpc-test-'));
  t.after(() => fs.rmSync(dir, { recursive: true, force: true }));
  const file = path.join(dir, 'test.pdf');
  fs.writeFileSync(file, bytes);
  return file;
}
function api(assets = [asset]) {
  const calls = [];
  const repos = {
    async getLatestRelease() { return { data: { id: 1, html_url: 'https://github.com/OI-wiki/OI-Wiki-export/releases/tag/test' } }; },
    async getReleaseByTag() { return { data: { id: 2, draft: false } }; },
    async createRelease() { calls.push('create'); return { data: { id: 2, draft: true } }; },
    async listReleaseAssets() {},
    async uploadReleaseAsset(options) { calls.push('upload'); return { data: { ...asset, id: 8, name: options.name } }; },
    async deleteReleaseAsset() { calls.push('delete'); },
    async updateReleaseAsset() { calls.push('rename'); },
    async updateRelease() { calls.push('update'); }
  };
  return { calls, repos, github: { rest: { repos },
    async paginate(_method, options) { return options.release_id === 1 ? [asset] : assets; }
  } };
}

test('verify accepts correct PDF and digest', async t => {
  assert.equal(await oi.verify(fixture(t), data.length, digest), digest);
});
test('verify rejects truncated downloads', async t => {
  await assert.rejects(oi.verify(fixture(t), data.length + 1), /size/);
});
test('verify rejects HTTP/HTML error payload', async t => {
  await assert.rejects(oi.verify(fixture(t, Buffer.from('error')), 5), /not a PDF/);
});
test('verify rejects incorrect hash', async t => {
  await assert.rejects(oi.verify(fixture(t), data.length, 'sha256:wrong'), /mismatch/);
});
test('matching metadata avoids unchanged daily downloads', async () => {
  const { github } = api();
  assert.equal((await oi.inspect({ github, context, core })).changed, 'false');
});
test('draft mirror is not considered published', async () => {
  const { github, repos } = api();
  repos.getReleaseByTag = async () => ({ data: { id: 2, draft: true } });
  assert.equal((await oi.inspect({ github, context, core })).changed, 'true');
});
test('permission errors are not swallowed as missing release', async () => {
  const { github, repos } = api();
  repos.getReleaseByTag = async () => { throw Object.assign(new Error('Forbidden'), { status: 403 }); };
  await assert.rejects(oi.inspect({ github, context, core }), /Forbidden/);
});
test('publish verifies upload before deleting previous asset', async t => {
  const { github, calls } = api([{ ...asset, digest: 'old' }]);
  await oi.publish({ github, context, core }, { file: fixture(t), digest, source: 'upstream' });
  assert.deepEqual(calls, ['upload', 'delete', 'rename', 'update']);
});
test('upload failure preserves existing public asset', async t => {
  const { github, calls, repos } = api([{ ...asset, digest: 'old' }]);
  repos.uploadReleaseAsset = async () => { calls.push('upload'); throw new Error('HTTP 500'); };
  await assert.rejects(oi.publish({ github, context, core }, { file: fixture(t), digest }), /HTTP 500/);
  assert.deepEqual(calls, ['upload']);
});
test('unchanged publish is idempotent', async t => {
  const { github, calls } = api();
  await oi.publish({ github, context, core }, { file: fixture(t), digest, source: 'upstream' });
  assert.deepEqual(calls, ['update']);
});
test('verified pending asset can resume without a second upload', async t => {
  const { github, calls } = api([{ ...asset, digest: 'old' },
    { ...asset, id: 8, name: `OI-wiki-${digest.slice(7)}.pdf` }]);
  await oi.publish({ github, context, core }, { file: fixture(t), digest, source: 'upstream' });
  assert.deepEqual(calls, ['delete', 'rename', 'update']);
});
test('PRs cannot call either publisher', async () => {
  const pr = { ...context, eventName: 'pull_request' };
  await assert.rejects(oi.publish({ github: {}, context: pr, core }), /default branch/);
  await assert.rejects(headers.publish({ github: {}, context: pr, core }), /default branch/);
});

function commitApi({ stale = false, identical = false, race = false, denied = false } = {}) {
  const calls = [];
  let advanced = false;
  const github = {
    rest: { git: {
      async getRef() { return { data: { object: { sha: advanced ? 'next' : 'head' } } }; },
      async getTree({ tree_sha }) { return { data: { truncated: false, tree: [
        { path: 'Header.hpp', sha: stale && tree_sha !== 'built' ? 'new-source' : 'source' },
        { path: 'headers.pdf', sha: identical ? headers.blobSha(data) : 'old-pdf' }
      ] } }; }
    } },
    async graphql(_query, { input }) {
      calls.push(input);
      if (denied) throw new Error('Permission denied');
      if (race && !advanced) { advanced = true; throw new Error('Head moved'); }
      return { createCommitOnBranch: { commit: { oid: 'new' } } };
    }
  };
  return { github, calls };
}
test('changed sources cannot be overwritten by a stale PDF build', async t => {
  const { github, calls } = commitApi({ stale: true });
  assert.equal(await headers.publish({ github, context, core }, fixture(t)), 'stale');
  assert.equal(calls.length, 0);
});
test('byte-identical headers PDF does not create a commit', async t => {
  const { github, calls } = commitApi({ identical: true });
  assert.equal(await headers.publish({ github, context, core }, fixture(t)), 'unchanged');
  assert.equal(calls.length, 0);
});
test('concurrent unrelated commits are retried with an updated expected head', async t => {
  const { github, calls } = commitApi({ race: true });
  assert.equal(await headers.publish({ github, context, core }, fixture(t)), 'published');
  assert.deepEqual(calls.map(c => c.expectedHeadOid), ['head', 'next']);
  assert.equal(calls[1].fileChanges.additions[0].path, 'headers.pdf');
});
test('real publication failures stay red', async t => {
  const { github, calls } = commitApi({ denied: true });
  await assert.rejects(headers.publish({ github, context, core }, fixture(t)), /Permission denied/);
  assert.equal(calls.length, 1);
});
