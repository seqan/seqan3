/* Enrico Seiler (enrico.seiler@fu-berlin.de), 2021 */

$(document).ready(function()
{
    var add_mlabel = function()
    {
        var memname_html = $(this).find('table.memname').prop('outerHTML');
        var new_html = "<table class=\"mlabels\"> <tbody> <tr> <td class=\"mlabels-left\"> " +
                        memname_html +
                        " </td> <td class=\"mlabels-right\"> <span class=\"mlabels\"\> " +
                        "</td> </tr> </tbody> </table>";

        $(this).find('div.memproto').html(new_html);
    }

    var add_labels = function(api_type, api_color, entry_obj, entry_type)
    {
        if (entry_type == "member")
        {
            $(this).find('span.mlabels').prepend('<span class="mlabel ' + api_type + '">' + api_type + '</span>');

            // get header above the item, append the api_type tag to the header and change the bg color
            prev = entry_obj;
            do prev = prev.previousSibling; while(prev && prev.nodeType !== 1);

            prev.innerContent = prev.innerContent + '  ' + '[' + api_type + '] ';

            if (api_type != "stable-api") // for now, leave stable api background in blue.
            {
                prev.style.backgroundImage = "none";
                prev.style.backgroundColor = api_color;
            }
        }
        else if (entry_type == "header")
        {
            $('h2.groupheader:contains("Detailed Description")').append('<span class="mlabel ' + api_type + ' header">' + api_type + '</span>');
        }
    }

    var select_and_add_label = function(entry_obj, entry_type)
    {
        var is_no_api = $(entry_obj).find('dl.no-api').length;
        var is_experimental_api = $(entry_obj).find('dl.experimental-api').length;
        var is_stable_api = $(entry_obj).find('dl.stable-api').length;
        var has_mlabels = $(this).find('span.mlabels').length;

        if (!has_mlabels)
            add_mlabel.call(this);

        if (is_no_api)
            add_labels.call(this, "no-api", '#ffb3b3', entry_obj, entry_type);
        else if (is_experimental_api)
            add_labels.call(this, "experimental-api", '#FBF0E0', entry_obj, entry_type);
        else if (is_stable_api)
            add_labels.call(this, "stable-api", '#e8ffe8', entry_obj, entry_type);
        else
            add_labels.call(this, "no-api", '#ffb3b3', entry_obj, entry_type);
    }

    // Badges for members
    $('div.memitem').each(function(index)
    {
        select_and_add_label.call(this, this, "member");
    });

    // Badge for "Detailed Description"
    select_and_add_label($('div.textblock'), "header");

});
