/// \file
/// \brief Robust estimation in presence of outliers with RANSAC and variants.
/// \details Two functions are proposed: Ransac (fix threshold for outliers and
/// minimize number of outliers) and LMedS (define outliers as quantile of worst
/// residuals and minimize the threshold). Template parameters should statisfy:
/// \li \c SampleSize number of data samples necessary to evaluate parameters.
/// \li \c DataInputIterator a forward iterator type over data.
/// \li \c Estimator functor of arguments
/// (FArray<DataInputIterator,SampleSize>,OutputIterator) putting in
/// \c OutputIterator the sets of parameters evaluated from samples in array
/// argument. OutputIterator::value_type must be \c Parameters.
/// \li \c Residual functor of arguments
/// (Parameters,DataInputIterator::value_type) returning a \c double error of
/// sample relative to the parameters.
/// \li \c Parameters set of parameters to estimate.
///
/// Normally, only the \c SampleSize template parameter must be precised, the
/// others are implicitely deduced from the arguments by the compiler.
#include <vector>

namespace Imagine {

	/// \addtogroup Optim
	/// @{

	/// \brief Get a random set of \a SampleSize samples out of \a n elements.
	/// \details \a out iterates over collection of type \a InputIterator.
	/// Only one pass is done through \a first, so it remains efficient in
	/// case of a forward iterator. If \a bUnique is set and the samples are not
	/// unique, return \a false. This is unlikely if SampleSize<<n, but in case the
	/// best is to try again.
	inline size_t size_tRandom(size_t a)
	{
		return (size_t)((a+.999)*doubleRandom());
	}
	template <int SampleSize, class InputIterator, class OutputIterator>
	bool getSample(InputIterator first, size_t n, OutputIterator out,
		bool bUnique=false)
	{
		FArray<std::pair<size_t,int>,SampleSize> index;
		for(int i=0; i < SampleSize; i++)
			index[i] = std::make_pair(size_tRandom(n-1), i);
		std::sort(index.begin(), index.end()); // Lexicographic order
		size_t prev=0;
		for(int i=0; i < SampleSize; i++) {
			if(bUnique && i != 0 && index[i].first == prev)
				return false;
			std::advance(first, index[i].first-prev);
			prev = index[i].first;
			OutputIterator o=out;
			std::advance(o, index[i].second);
			*o = first;
		}
		return true;
	}

	/// Number of samplings to have at least \a minProba probability of absence of
	/// outlier in a sample of \a SampleSize elements.
	inline int getNumSamples(double minProba, double outlierRatio, int SampleSize)
	{
		return int( std::log(1.-minProba) /
			std::log(1.-std::pow(1.-outlierRatio, SampleSize)) );
	}

	/// \brief The famous Random Sample Consensus algorithm (Fischler&Bolles 1981).
	/// \details The number of tests is reevaluated down as soon as more inliers are
	/// found. Return the number of inliers.
	template <int SampleSize,  class DataInputIterator, class Estimator,
	class Residual, class Parameters>
		std::size_t Ransac(DataInputIterator first, DataInputIterator last,
		const Estimator& estimator, const Residual& residual,
		Parameters& bestParameters,
		double outlierThreshold,
		double maxOutlierRatio=0.5,
		double minProba=0.99)
	{
		std::size_t bestInliers = 0;
		std::size_t n = std::distance(first,last);

		FArray<DataInputIterator,SampleSize> sample;

		// Required number of iterations is evaluated from outliers ratio
		int N = (SampleSize<=n)?
			getNumSamples(minProba, maxOutlierRatio, SampleSize): 0;
		for (int i=0; i < N; i++) {
			// Try until getting a set of SampleSize unique samples
			while(! getSample<SampleSize>(first, n, sample.begin(), true)) {}
			// Estimate parameters: the solutions are stored in a vector
			std::vector<Parameters> paramVec;
			estimator(sample, std::back_inserter(paramVec));

			// Now test the solutions on the whole data
			typename std::vector<Parameters>::iterator it = paramVec.begin();
			for(; it != paramVec.end(); ++it) {
				std::size_t inliers = 0;
				for(DataInputIterator it2=first; it2 != last; ++it2)
					if(residual(*it,*it2) < outlierThreshold)
						++ inliers;

				// Store best solution and reestimate number of tests
				if (inliers > bestInliers) {
					bestInliers = inliers;
					bestParameters = *it;
					double outlierRatio = (n - inliers)/(double)n;
					if(outlierRatio < maxOutlierRatio)
						N = getNumSamples(minProba, outlierRatio, SampleSize);
				}
			}
		}
		return bestInliers;
	}

	///Modified Ransac: return solution with less residual
	template <int SampleSize,  class DataInputIterator, class Estimator,
	class Residual, class Parameters>
		std::size_t ModifiedRansac(DataInputIterator first, DataInputIterator last,
		const Estimator& estimator, const Residual& residual,
		Parameters& bestParameters,
		double outlierThreshold,
		double maxOutlierRatio=0.5,
		double minProba=0.99)
	{
		std::size_t bestInliers = 0;
		std::size_t n = std::distance(first,last);

		FArray<DataInputIterator,SampleSize> sample;

		// Required number of iterations is evaluated from outliers ratio
		int N = (SampleSize<=n)?
			getNumSamples(minProba, maxOutlierRatio, SampleSize): 0;
		for (int i=0; i < N; i++) {
			// Try until getting a set of SampleSize unique samples
			while(! getSample<SampleSize>(first, n, sample.begin(), true)) {}
			// Estimate parameters: the solutions are stored in a vector
			std::vector<Parameters> paramVec;
			estimator(sample, std::back_inserter(paramVec));

			// Now test the solutions on the whole data
			double min_sumresidual = std::numeric_limits<double>::infinity();
			typename std::vector<Parameters>::iterator it = paramVec.begin();
			for(; it != paramVec.end(); ++it) {
				std::size_t inliers = 0;
				double sum_residual = 0;
				double res;
				for(DataInputIterator it2=first; it2 != last; ++it2)
				{
    	 		res =	residual(*it,*it2);
					if(res < outlierThreshold)
					{
						++ inliers;
						sum_residual += res;
				  }
        }
        
				// Store best solution and reestimate number of tests
				if (inliers > bestInliers || (inliers == bestInliers && sum_residual < min_sumresidual)) {
					bestInliers = inliers;
					min_sumresidual = sum_residual;
					bestParameters = *it;
					double outlierRatio = (n - inliers)/(double)n;
					if(outlierRatio < maxOutlierRatio)
						N = getNumSamples(minProba, outlierRatio, SampleSize);
				}
			}
		}
		return bestInliers;
	}



	/// \brief Variant of RANSAC using Least Median of Squares.
	/// \details Instead of using a fixed threshold to distinguish inlier/outlier,
	/// find the threshold at 1-\a outlierRatio quantile of residuals and keep the
	/// parameters minimizing this threshold. The final threshold
	/// returned in \a outlierThreshold is a multiple of this and can be used to
	/// filter out outliers.
	template <int SampleSize, class DataInputIterator,
	class Estimator, class Residual, class Parameters>
		double LMedS(DataInputIterator first, DataInputIterator last,
		const Estimator& estimator, const Residual& residual,
		Parameters& bestParameters,
		double* outlierThreshold=NULL,
		double outlierRatio=0.5,
		double minProba=0.99)
	{
		std::size_t n = std::distance(first,last); // Data size
		std::vector<double> residuals(n); // Array for storing residuals

		FArray<DataInputIterator,SampleSize> sample;
		double bestMedian = std::numeric_limits<double>::max();

		// Required number of iterations is evaluated from outliers ratio
		const int N = (SampleSize<n)?
			getNumSamples(minProba, outlierRatio, SampleSize): 0;
		for (size_t i=0; i < N; i++) {
			// Try until getting a set of SampleSize unique samples
			while(! getSample<SampleSize>(first, n, sample.begin(), true)) {}
			// Estimate parameters: the solutions are stored in a vector
			std::vector<Parameters> paramVec;
			estimator(sample, std::back_inserter(paramVec));

			// Now test the solutions on the whole data
			typename std::vector<Parameters>::iterator it = paramVec.begin();
			for(; it != paramVec.end(); ++it) {
				// Compute residuals
				std::vector<double>::iterator out = residuals.begin();
				for(DataInputIterator it2 = first; it2 != last; ++it2)
					*out++ = residual(*it,*it2);

				// Compute median
				std::vector<double>::iterator itMedian = residuals.begin() +
					std::size_t( n*(1.-outlierRatio) );
				std::nth_element(residuals.begin(), itMedian, residuals.end());
				double median = *itMedian;

				// Store best solution
				if(median < bestMedian) {
					bestMedian = median;
					bestParameters = *it;
				}
			}
		}

		/* this array of precomputed values corresponds to the inverse
		cumulative function for a normal distribution. For more information
		consult the litterature (Robust Regression for Outlier Detection,
		rouseeuw-leroy). The values are computed for each 5% */
		static const double ICDF[21] = {
			1.4e16, 15.94723940, 7.957896558, 5.287692054, 
			3.947153876, 3.138344200, 2.595242369, 2.203797543, 
			1.906939402, 1.672911853, 1.482602218, 1.323775627, 
			1.188182950, 1.069988721, 0.9648473415, 0.8693011162, 
			0.7803041458, 0.6946704675, 0.6079568319,0.5102134568,
			0.3236002672
		};

		// Evaluate the outlier threshold
		if(outlierThreshold) {
			double sigma = ICDF[int((1.-outlierRatio)*20.)] *
				(1. + 5. / double(n - SampleSize));
			*outlierThreshold = (double)(sigma * sigma * bestMedian * 4.);
		}

		return bestMedian;
	}

	/// \brief Keep only inliers in data set \a first to \a last.
	/// \details The returned iterator marks the new end of the list. For example,
	/// if data is in a \code std::list ls \endcode, you can call
	/// \code ls.erase(FilterOutliers(ls.begin(),ls.end(),res,param,th), ls.end());
	/// \endcode
	template <class DataIterator, class Residual, class Parameters>
	DataIterator FilterOutliers(DataIterator first, DataIterator last,
		Residual residual,
		Parameters param, double outlierThreshold)
	{
		DataIterator it=first;
		for(; first != last; first++)
			if(residual(param,*first) < outlierThreshold)
				*it++ = *first;
		return it;
	}
	///@}

}
