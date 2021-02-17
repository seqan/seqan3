/* Enrico Seiler (enrico.seiler@fu-berlin.de), 2021 */

$(document).ready(function()
{
    $('div.memitem').each(function(index)
    {
        var is_no_api = $(this).find('dl.no-api').length;
        var is_experimental_api = $(this).find('dl.experimental-api').length;
        var is_stable_api = $(this).find('dl.stable-api').length;
        var has_mlabels = $(this).find('span.mlabels').length;

        if (!has_mlabels)
        {
            var memname_html = $(this).find('table.memname').prop('outerHTML');
            var new_html = "<table class=\"mlabels\"> <tbody> <tr> <td class=\"mlabels-left\"> " +
                           memname_html +
                           " </td> <td class=\"mlabels-right\"> <span class=\"mlabels\"\> " +
                           "</td> </tr> </tbody> </table>";
            $(this).find('div.memproto').html(new_html);
        }

        if (is_no_api)
            $(this).find('span.mlabels').prepend('<span class="mlabel no-api">no-api</span>');
        else if (is_experimental_api)
            $(this).find('span.mlabels').prepend('<span class="mlabel experimental-api">experimental-api</span>');
        else if (is_stable_api)
            $(this).find('span.mlabels').prepend('<span class="mlabel stable-api">stable-api</span>');
        else
            $(this).find('span.mlabels').prepend('<span class="mlabel no-api">no-api</span>');
    });
});
