# Contributing {#about_contributing}

<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0 -->

First of all, thanks for wanting to contribute to SeqAn! Community is important to us and we strive to maintain a great
culture and atmosphere. Please have a look at our [Code of Conduct](\ref about_code_of_conduct).

The following is a guide that helps stream-line the process of adding changes. If you haven't contributed to SeqAn
before, don't worry about getting something wrong, that's absolutely OK. However, following this guide closely
will reduce work for all of us and increases the chance of your changes being merged quickly.

# Pull request workflow

The overall workflow of contributing changes to this repository is:

  1. Create an issue with bug report or feature request.
  2. Wait for comments from the team and then collectively decide on the best way to solve the issue.
  3. Fork the repository and create a branch on your fork.
  4. Commit the solution to your local branch.
  5. Create a pull-request to our main branch from your branch.
  6. Wait for reviews and resolve reviewer comments.
  7. Wait for the PR to be merged.

## Creating an issue

Go to https://github.com/seqan/seqan3/issues and check if someone else has already reported the same issue. If not,
create a new issue and follow the template.

We usually reply quickly to issues; if there is no reply within a week (and it's not prime vacation season in Europe)
feel free to reply to the issue and explicitly mention `@seqan/core`.

## Creating commits

If you have never used the Fork & Pull Model, have a look at [this guide](https://guides.github.com/introduction/flow/).

When creating commits to your branch and before requesting a pull, please ensure the following:

  * Proper **naming** of commits (see below).
  * Content other than the first line is irrelevant in the commit message.
  * The number of commits should be minimal and commits shouldn't revert/change previous ones.
    * You can use separate commits for implementation/tests/documentation or for individual features/fixes.
    * Do not add additional commits a la "Resolve reviewer comments", instead
      [rewrite your old commits](https://thoughtbot.com/blog/git-interactive-rebase-squash-amend-rewriting-history).
  * When you add a feature, add a note to the CHANGELOG.md as part of your changes.

Commit message prefixes:
```
[FEATURE] Whenever you implement something new and shiny
[FIX]     Whenever you fix some wrong code in the source
[DOC]     Whenever you do something only(!) related to the documentation
[INFRA]   Whenever you change something of the build system or CI related
[TEST]    Whenever you do something related to the tests (unit or benchmark)
[MISC]    Whenever it does not fit to any of the above
```

E.g. `[FEATURE] add support for FOO` or `[FIX] divide by zero in alignment`.

**Alternatively, if your changes are so small that they would fit into one commit, don't worry about any of the above,
just state "Please squash commits" in the PR and we will take care of the rest.**

## Opening a pull-request ("PR")

When you are done with committing changes to your branch and you have tested your changes, open a pull request.
We have continuous integration in place that should inform you of test failures. Please try to resolve any
breakage that your pull request introduces.

[Here is a guide for setting up unit tests locally.](https://docs.seqan.de/seqan3/main_user/setup_tests.html)

If there are test failures that you don't understand, clearly indicate that you have seen the errors, but cannot resolve
them – then the first reviewer will have a look at them. Otherwise the PR will be treated as still being
work-in-progress.

## The review process

### Request first review

After opening a PR it goes through the review process. This is a two-step process, first a regular member of the team
needs to approve your changes, then one of the project owners (`@h-2` or `@rrahn`) needs to approve and merge it. If you
are a collaborator of the project and know the SeqAn team member best suitable for the review you can request a review
yourself. Otherwise wait for a reviewer to be assigned to the PR. Do not request reviews from a project owner at this
point. Always only request review from one person at a time.

### Resolve comments

If you have not contributed to SeqAn before, you will likely receive a lot of comments on the style of your code. This
is not a sign that we don't appreciate your contribution, but we have high code quality standards and everything
must conform to our [style guide](https://github.com/seqan/seqan3/wiki#library-coding-guide). In the future we hope
to automate this step.
After you have received the first review, go through all comments and clarify any suggestions and requests. When
you are confident that you know how to address all comments, implement the changes and update your PR.

### Re-request review

Now re-request a review from the previous reviewer either by clicking "re-request review" or by writing "please review
again" as a comment. **This step is important**, otherwise your PR is treated as still being work-in-progress! The
reviewer will go through their first review and mark all previous comments as resolved if they think you resolved them
(**please do not mark another person's comments as resolved**). They will then add a second review or approve your
changes. Repeat the last steps until you get approval.

### Request project owner review

Once your changes are approved by one reviewer, request a second review from a project owner (or wait to be assigned).
They will go through everything briefly again and make sure that nothing escaped the first reviewer's notice.
Ultimately, they will approve your changes and merge your branch.

### Post scriptum

We try to review all PRs as soon as possible. If you don't get a response within a week, ping the requested reviewer
politely. Working with git and understanding the review process can be challenging at first. Please don't interpret our
insistence on certain changes as a rejection of your contribution – we try to do our best in maintaining high quality
software and extensive reviews are a part of this process. If at any time you feel unable to resolve a comment, just
ask the reviewer for guidance. We can also push changes to your branch to fix things, but we will only do so if you
ask us to.
