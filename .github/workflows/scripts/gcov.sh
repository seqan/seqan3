#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0


args=${@/--branch-counts/""}
args=${args/--branch-probabilities/""}

exec gcov $args
