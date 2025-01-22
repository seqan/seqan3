/*  SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
    SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
    SPDX-License-Identifier: BSD-3-Clause
*/
/* Jongkyu Kim (j.kim@fu-berlin.de), 2016.01.12
   Adaptations by Enrico Seiler (enrico.seiler@fu-berlin.de), 2020 */

function changeVersion(form_id)
{
    // Get the base url without version information, e.g. "https://docs.seqan.de/seqan"
    var current_script_url = document.scripts[document.scripts.length - 1].src;
    var base_url = current_script_url.split('/').slice(0, -2).join('/');

    // Get the current page, e.g. "index.html"
    var full_url = window.top.location.href;
    var current_page = full_url.substring(full_url.lastIndexOf("/") + 1);

    // Get the selected version
    var form = document.getElementById(form_id);
    var version = form.options[form.selectedIndex].value;

    // Check if the current page is valid with the selected version
    var proposed_url = base_url + '/' + version + '/' + current_page;
    var request = new XMLHttpRequest();
    request.open('GET', proposed_url, false);
    request.send();
    // If the URL is invalid, redirect to main page of the selected version
    // If htaccess is configured to redirect invalid URLs to the base domain,
    // no 404 is returned, hence the second condition
    if (request.status === 404 || request.responseURL == window.location.origin + '/')
    {
        proposed_url = base_url + '/' + version;
    }

    // Load the proper page
    window.top.location.href = proposed_url;
}

function addVersionSelection(arr)
{
    // add HTMLs
    var version_select = document.createElement("select");
    version_select.setAttribute("id","version_select");
    version_select.setAttribute("style","color: var(--page-foreground-color);");
    document.getElementById("list_bottom_right").appendChild(version_select);

    version_select.addEventListener("change", function(){changeVersion(this.id);}, false);

    // current selection is..
    cur_sel = window.location.pathname.split("/").at(-2);

    for(i=0; i < arr.length; ++i)
    {
        var op = document.createElement("option");
        op.value = arr[i];
        op.text = arr[i];
        op.selected = ( arr[i] == cur_sel ) ? true : false;
        version_select.add(op);
    }
}

// get JSON data & add selection form
var request = new XMLHttpRequest();
request.open("GET", "version.php", true);
request.setRequestHeader("Content-type", "application/json");
request.onreadystatechange = function()
{
    if( request.readyState == 4 && request.status == 200 )
    {
        var response = JSON.parse(request.responseText);
        addVersionSelection(response); // add selection form
    }
}
request.send();
