<!-- SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

How to setup vercel.com:

* In vercel.com go create a new project.
* Import it from the git repository.
* Select personal account
* Give it a Project Name and select Framework Preset to Other
* Set Output Directory to "./export"
* Set Build Command to "./test/documentation/.vercel/build.sh"
* This option needs to be changed in the settings of the project:
  * Set Install Command to "./test/documentation/.vercel/install.sh"
