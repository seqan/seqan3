<?php
/* Jongkyu Kim(j.kim@fu-berlin.de), 2016.01.12
   Adaptations by Enrico Seiler (enrico.seiler@fu-berlin.de), 2020 */
$LOCALDIR = "../";

$files = scandir($LOCALDIR, SCANDIR_SORT_DESCENDING);
$list = array();
foreach( $files as $file  )
{
    if( strpos($file, "master_") !== FALSE ) // Skip directories starting with "master_"
        continue;
    if( $file == "master") // Skip SeqAn2 "master" documentation
        continue;
    if( $file == "index.html") // Skip seqan/index.html
        continue;
    if( strpos($file, "develop") !== FALSE )  // Skip SeqAn2 "develop" and "develop_" documentation
        continue;
    if( $file[0] == "1") // Skip SeqAn1 versioned documentation
        continue;
    if( $file[0] == "2") // Skip SeqAn2 versioned documentation
        continue;
    if( strpos($file, "learning-resources") !== FALSE ) // Skip learning resources
        continue;
    if( $file[0] == "." ) // Skip current directory
        continue;

    array_push($list, $file);
}

header("Content-Type: application/json");
echo json_encode($list, JSON_PRETTY_PRINT);
?>
