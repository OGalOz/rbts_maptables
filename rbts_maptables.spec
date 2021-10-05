/*
A KBase module: rbts_maptables
*/

module rbts_maptables {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_rbts_maptables(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

};
