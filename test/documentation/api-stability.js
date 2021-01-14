/* Enrico Seiler (enrico.seiler@fu-berlin.de), 2021 */

$(document).ready(function()
{
    $('div.memitem').each(function(index) 
    {
        var is_no_api = $(this).find('dl.no-api').length;
        var is_experimental_api = $(this).find('dl.experimental-api').length;
        var is_stable_api = $(this).find('dl.stable-api').length;

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
