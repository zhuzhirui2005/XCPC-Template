'use strict';
const fs = require('node:fs');
const path = require('node:path');
const { createHash } = require('node:crypto');
const TAG = 'oi-wiki-latest';
const NAME = 'OI-wiki.pdf';
const matches = (asset, size, digest) => Boolean(asset && digest &&
  asset.state === 'uploaded' && asset.size === size && (asset.digest || asset.label) === digest);

async function inspect({ github, context, core }) {
  const upstream = { owner: 'OI-wiki', repo: 'OI-Wiki-export' };
  const { data: release } = await github.rest.repos.getLatestRelease(upstream);
  const assets = await github.paginate(github.rest.repos.listReleaseAssets, {
    ...upstream, release_id: release.id, per_page: 100
  });
  const asset = assets.find(a => a.name === NAME);
  if (!asset || asset.state !== 'uploaded' || !Number.isSafeInteger(asset.size) || asset.size <= 0) {
    throw new Error('Cannot find a complete OI-wiki.pdf in the latest release');
  }
  const digest = asset.digest || '';
  if (digest && !/^sha256:[a-f0-9]{64}$/.test(digest)) throw new Error('Unsupported upstream digest');
  let changed = true;
  try {
    const { data: mirror } = await github.rest.repos.getReleaseByTag({ ...context.repo, tag: TAG });
    const mirrored = await github.paginate(github.rest.repos.listReleaseAssets, {
      ...context.repo, release_id: mirror.id, per_page: 100
    });
    changed = mirror.draft || !matches(mirrored.find(a => a.name === NAME), asset.size, digest);
  } catch (error) {
    if (error.status !== 404) throw error;
  }
  const outputs = { url: asset.browser_download_url, size: String(asset.size), digest,
    source_url: release.html_url, changed: String(changed) };
  for (const [key, value] of Object.entries(outputs)) core.setOutput(key, value);
  core.info(changed ? 'New or unverified upstream PDF.' : 'Already current; no PDF download needed.');
  return outputs;
}

async function verify(file, size, expected = '') {
  if (!Number.isSafeInteger(Number(size)) || Number(size) <= 0 || fs.statSync(file).size !== Number(size)) {
    throw new Error('Downloaded PDF size does not match the release asset');
  }
  const handle = fs.openSync(file, 'r');
  const header = Buffer.alloc(5);
  try { fs.readSync(handle, header, 0, 5, 0); } finally { fs.closeSync(handle); }
  if (header.toString('ascii') !== '%PDF-') throw new Error('Downloaded file is not a PDF');
  const hash = createHash('sha256');
  for await (const chunk of fs.createReadStream(file)) hash.update(chunk);
  const digest = `sha256:${hash.digest('hex')}`;
  if (expected && expected !== digest) throw new Error('Downloaded PDF SHA-256 mismatch');
  return digest;
}

async function publish({ github, context, core }, options = {}) {
  const branch = context.payload.repository.default_branch;
  if (context.eventName === 'pull_request' || context.ref !== `refs/heads/${branch}`) {
    throw new Error('Release publication is only allowed on the default branch');
  }
  const file = options.file || path.join(process.env.RUNNER_TEMP, NAME);
  const digest = options.digest || process.env.PDF_DIGEST;
  const source = options.source || process.env.SOURCE_RELEASE_URL;
  if (!/^sha256:[a-f0-9]{64}$/.test(digest || '')) throw new Error('Missing verified PDF digest');
  const size = fs.statSync(file).size;
  await verify(file, size, digest); // Also verify the artifact after its cross-job transfer.
  const data = fs.readFileSync(file);
  const repo = context.repo;
  const body = [
    `Automatically mirrored from [OI-Wiki-export](${source}).`, '',
    'Download `OI-wiki.pdf` below. This is a release asset, not a Git blob.', '',
    `SHA-256: \`${digest.slice(7)}\``
  ].join('\n');
  let release;
  try {
    ({ data: release } = await github.rest.repos.getReleaseByTag({ ...repo, tag: TAG }));
  } catch (error) {
    if (error.status !== 404) throw error;
    ({ data: release } = await github.rest.repos.createRelease({
      ...repo, tag_name: TAG, target_commitish: context.sha,
      name: 'OI-Wiki PDF (latest)', body, draft: true, make_latest: 'false'
    }));
  }
  const assets = await github.paginate(github.rest.repos.listReleaseAssets, {
    ...repo, release_id: release.id, per_page: 100
  });
  const current = assets.find(a => a.name === NAME);
  if (!matches(current, size, digest)) {
    // Complete and verify the upload before replacing the public download.
    const pendingName = `OI-wiki-${digest.slice(7)}.pdf`;
    let pending = assets.find(a => a.name === pendingName);
    if (!matches(pending, size, digest)) {
      if (pending) await github.rest.repos.deleteReleaseAsset({ ...repo, asset_id: pending.id });
      ({ data: pending } = await github.rest.repos.uploadReleaseAsset({
        ...repo, release_id: release.id, name: pendingName, label: digest,
        headers: { 'content-type': 'application/pdf', 'content-length': size }, data
      }));
      if (!matches(pending, size, digest)) throw new Error('Uploaded PDF asset failed verification');
    }
    if (current) await github.rest.repos.deleteReleaseAsset({ ...repo, asset_id: current.id });
    await github.rest.repos.updateReleaseAsset({ ...repo, asset_id: pending.id, name: NAME, label: digest });
  }
  if (release.draft || release.body !== body) {
    await github.rest.repos.updateRelease({ ...repo, release_id: release.id, body, draft: false, make_latest: 'false' });
  }
  core.info(`Published ${NAME}: ${digest}`);
}
module.exports = { inspect, verify, publish, matches };
