#include "iset.h"

SEXP get_box_bounds (SEXP runs, SEXP reflevels, SEXP adex, SEXP chrs, SEXP starts, SEXP ends) {
    BEGIN_RCPP

    // Setting up input objects.
    const Rcpp::IntegerVector Runs(runs), Adex(adex);
    const size_t nlevels=Runs.size();
    const size_t npts=std::accumulate(Runs.begin(), Runs.end(), 0);
    if (npts!=Adex.size()) { 
        throw std::runtime_error("anchor index length should be equal to sum of factor runs");
    }
    
    const Rcpp::IntegerVector Chrs(chrs), Starts(starts), Ends(ends);
    const int nregs=Chrs.size();
    if (nregs!=Starts.size() || nregs!=Ends.size()) {
        throw std::runtime_error("chromosome/start/end vectors should have the same length"); 
    }

    for (const auto& A : Adex) {
        if (A < 0 || A >= nregs) {
            throw std::runtime_error("anchor index out of range of 'regions(x)'");
        }
    }

    // Setting up output constructs.
    Rcpp::IntegerVector out_chr(nlevels), out_start(nlevels), out_end(nlevels);
    Rcpp::IntegerVector::const_iterator aIt=Adex.begin();
    Rcpp::IntegerVector::iterator ocIt=out_chr.begin(), osIt=out_start.begin(), oeIt=out_end.begin();
    int run_counter=0;

    for (const auto& rlen : Runs) {
        const int& curA=*aIt;
        int& curchr=(*ocIt = Chrs[curA]);
        int& curstart=(*osIt = Starts[curA]);
        int& curend=(*oeIt = Ends[curA]);
        ++aIt;        

        for (int i=1; i<rlen; ++i) {
            const int& nextA=*aIt;

            if (Chrs[nextA]!=curchr) {
                Rcpp::StringVector mylevels(reflevels); // Coerces it to a string, if it wasn't already.
                if (mylevels.size() <= run_counter) {
                    throw std::runtime_error("insufficient levels supplied for the given number of runs");
                }
                std::stringstream err_out;
                err_out << "multiple chromosomes for group '" << Rcpp::as<std::string>(mylevels[run_counter]) << "'";
                throw std::runtime_error(err_out.str());
            }
            if (Starts[nextA] < curstart) { 
                curstart = Starts[nextA];
            }
            if (Ends[nextA] > curend) { 
                curend = Ends[nextA];
            }
            ++aIt;
        }

        ++ocIt;
        ++osIt;
        ++oeIt;
        ++run_counter;
    }
   
    return Rcpp::List::create(out_chr, out_start, out_end);
    END_RCPP
}

