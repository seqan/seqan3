<?php
/* Jongkyu Kim(j.kim@fu-berlin.de), 2016.01.12 */
$LOCALDIR = "../";

$files = scandir($LOCALDIR, SCANDIR_SORT_DESCENDING);
$list = array();
foreach( $files as $file  )
{
    if( strpos($file, "master_") !== FALSE )
        continue;
    if( $file == "master")
        continue;
    if( strpos($file, "develop") !== FALSE )
        continue;
    if( $file[0] == "1")
        continue;
    if( $file[0] == "2")
        continue;
    if( strpos($file, "learning-resources") !== FALSE )
        continue;
    if( $file[0] == "." )
        continue;

    array_push($list, $file);
}

header("Content-Type: application/json");
echo json_encode($list, JSON_PRETTY_PRINT);
?>
