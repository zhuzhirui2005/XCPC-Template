'use strict';
const fs = require('node:fs');
const { createHash } = require('node:crypto');

const isInput = p => (!p.includes('/') && p.endsWith('.hpp')) || [
  '.github/scripts/headers.py', '.github/workflows/generate-headers-pdf.yml'
].includes(p);
const signature = tree => JSON.stringify(tree.filter(x => isInput(x.path))
  .map(x => [x.path, x.sha]).sort((a, b) => a[0].localeCompare(b[0], 'en')));
const blobSha = data => createHash('sha1').update(`blob ${data.length}\0`).update(data).digest('hex');

async function publish({ github, context, core }, file = process.env.PDF_PATH) {
  const branch = context.payload.repository.default_branch;
  if (context.eventName === 'pull_request' || context.ref !== `refs/heads/${branch}`) {
    throw new Error('PDF commits are only allowed on the default branch');
  }
  const data = fs.readFileSync(file);
  if (data.subarray(0, 5).toString() !== '%PDF-' || data.length >= 100 * 1024 * 1024) {
    throw new Error('Invalid or oversized headers.pdf');
  }
  const repo = context.repo;
  const getTree = async sha => {
    const { data: tree } = await github.rest.git.getTree({ ...repo, tree_sha: sha, recursive: '1' });
    if (tree.truncated) throw new Error('Cannot verify a truncated repository tree');
    return tree.tree;
  };
  const inputs = signature(await getTree(context.sha));
  for (let attempt = 0; attempt < 3; attempt++) {
    const { data: ref } = await github.rest.git.getRef({ ...repo, ref: `heads/${branch}` });
    const head = ref.object.sha;
    const tree = await getTree(head);
    if (signature(tree) !== inputs) {
      core.warning('PDF inputs changed during this build; not committing a stale PDF.');
      return 'stale';
    }
    if (tree.find(x => x.path === 'headers.pdf')?.sha === blobSha(data)) {
      core.info('headers.pdf is unchanged; skipping commit.');
      return 'unchanged';
    }
    try {
      await github.graphql(`mutation($input: CreateCommitOnBranchInput!) {
        createCommitOnBranch(input: $input) { commit { oid } }
      }`, { input: {
        branch: { repositoryNameWithOwner: `${repo.owner}/${repo.repo}`, branchName: branch },
        expectedHeadOid: head,
        message: { headline: 'chore: update headers PDF' },
        fileChanges: { additions: [{ path: 'headers.pdf', contents: data.toString('base64') }] }
      } });
      core.info('Updated headers.pdf without rewriting branch history.');
      return 'published';
    } catch (error) {
      const { data: latest } = await github.rest.git.getRef({ ...repo, ref: `heads/${branch}` });
      if (latest.object.sha === head) throw error;
      // A concurrent commit won the compare-and-swap; recheck its inputs.
    }
  }
  throw new Error('Branch kept changing; refusing to overwrite a newer commit');
}
module.exports = { publish, signature, blobSha };
