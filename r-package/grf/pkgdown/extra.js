// From https://github.com/microsoft/LightGBM/blob/master/R-package/pkgdown/extra.js
$(function() {
    if(window.location.pathname.toLocaleLowerCase().indexOf('/reference') != -1) {
        /* Replace '/R/' with '/r-package/grf/R/' in all external links to .R files of grf GitHub repo */
        $('a[href^="https://github.com/grf-labs/grf/blob/master/R"][href*=".R"]').attr('href', (i, val) => { return val.replace('/R/', '/r-package/grf/R/'); });
    }
});
