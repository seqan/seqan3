<?php
/*
    SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
    SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
    SPDX-License-Identifier: BSD-3-Clause
*/
/* Jongkyu Kim(j.kim@fu-berlin.de), 2016.01.12
   Adaptations by Enrico Seiler (enrico.seiler@fu-berlin.de), 2020 */
$LOCALDIR = "../";

$files = scandir($LOCALDIR, SCANDIR_SORT_DESCENDING);
$list = array();
foreach( $files as $file  )
{
    if( !is_dir("$LOCALDIR/$file") ) // Skip files
        continue;
    if( strpos($file, "hidden_") !== FALSE ) // Skip directories starting with "hidden_"
        continue;
    if( $file[0] == "." ) // Skip current directory
        continue;

    array_push($list, $file);
}

header("Content-Type: application/json");
echo json_encode($list, JSON_PRETTY_PRINT);
?>
