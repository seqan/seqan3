## api.yml

Runs the [API-Stability](https://github.com/seqan/seqan3/blob/master/test/api_stability/README.md) test:
  * Every Sunday on master
  * Manually on any branch

In case of failure, creates an issue containing error logs.

## cancel.yml

Cancels ongoing and queued runs for a branch if a new push is made:
  * Any push to a branch
  * Any push to a pull request

This worklow needs to run on `pull_request_target` instead of `pull_request`.
This means that this workflow runs on the base branch (i.e., `seqan3/master`) which is necessary to have the write
permissions needed to cancel a run.

## ci_linux.yml

Runs the unit test CI on Ubuntu runners:
  * Any push to a master or release branch
  * Any push to a pull request

Also produces code coverage.

## ci_macos.yml

Runs the unit test CI on macOS runners:
  * Any push to a master or release branch
  * Any push to a pull request

In contrast to Ubuntu runners, there are 3 cores available.
Since macOS runs take longer and there generally less workers available, this workflow should be kept short.

## ci_misc.yml

Runs various CI on Ubuntu runners:
  * Any push to a master or release branch
  * Any push to a pull request

This workflow includes performance and snippet tests, CMake related CI, header test, and the documentation build.

## documentation.yml

Builds documentation and uploads it to https://docs.seqan.de/:
  * Any push to master
  * Manually on any branch

Uses repository secrets for the target server and authentication.

## latest_libraries.yml

Will update all submodules, gbenchmark and googletest to respective master versions and runs the test suite on all
compilers:
  * Every Sunday on master
  * Manually on any branch

In case of failure, creates an issue containing error logs.

## ram_usage.yml

Tracks the RAM-Usage when compiling and generates a CSV as artefact:
  * Manually on any branch

This workflow needs two inputs:
  * The GCC version (e.g., `11`)
  * The CXX_FLAGS (optional)
