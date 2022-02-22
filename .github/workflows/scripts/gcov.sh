#!/usr/bin/env bash

args=${@/--branch-counts/""}
args=${args/--branch-probabilities/""}

exec gcov $args
