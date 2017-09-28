#include "iset.h"

SEXP get_box_bounds (SEXP ids, SEXP reflevels, SEXP adex, SEXP chrs, SEXP starts, SEXP ends) {
    BEGIN_RCPP

    // Setting up input objects.
    const Rcpp::IntegerVector _ids(ids),_adex(adex);
    const int npts=_ids.size();
    if (npts!=_adex.size()) { 
        throw std::runtime_error("anchor index and grouping vectors are not of same length"); 
    }

    const Rcpp::StringVector _reflevels(reflevels);
    const int nlevels=_reflevels.size();
    
    const Rcpp::IntegerVector _chrs(chrs), _starts(starts), _ends(ends);
    const int nregs=_chrs.size();
    if (nregs!=_starts.size() || nregs!=_ends.size()) {
        throw std::runtime_error("chromosome/start/end vectors should have the same length"); 
    }

    int ngroups=(npts > 0 ? 1 : 0);
    for (int i=1; i<npts; ++i) {
        if (_ids[i]!=_ids[i-1]) { ++ngroups; }
    }

    // Setting up output constructs.
    Rcpp::IntegerVector out_fac(ngroups), out_chr(ngroups), out_start(ngroups), out_end(ngroups);
    if (ngroups > 0) {
        out_fac[0]=_ids[0];
        out_chr[0]=_chrs[_adex[0]];
        out_start[0]=_starts[_adex[0]];
        out_end[0]=_ends[_adex[0]];
    }

    int curgroup=0;
    for (int i=1; i<npts; ++i) {
        const int& curdex=_adex[i];
        if (curdex < 0 || curdex >= nregs) { 
            throw std::runtime_error("anchor index is out of bounds"); 
        }

        const int& curid=_ids[i];
        if (curid!=_ids[i-1]) {
            ++curgroup;
            if (curgroup >= ngroups) { 
                throw std::runtime_error("internal indexing error for groups"); 
            }
            out_fac[curgroup]=curid;
            out_chr[curgroup]=_chrs[curdex];
            out_start[curgroup]=_starts[curdex];
            out_end[curgroup]=_ends[curdex];
        } else {
            if (out_chr[curgroup]!=_chrs[curdex]) {
                if (curid >= nlevels) {
                    throw std::runtime_error("grouping index outside range of reference levels");
                }
                std::stringstream err_out;
                err_out << "multiple chromosomes for group '" << Rcpp::as<std::string>(_reflevels[curid]) << "'";
                throw std::runtime_error(err_out.str());
            }
            if (out_start[curgroup] > _starts[curdex]) { 
                out_start[curgroup] = _starts[curdex]; 
            }
            if (out_end[curgroup] < _ends[curdex]) { 
                out_end[curgroup] = _ends[curdex]; 
            }
        }
    }
    
    return Rcpp::List::create(out_fac, out_chr, out_start, out_end);
    END_RCPP
}

