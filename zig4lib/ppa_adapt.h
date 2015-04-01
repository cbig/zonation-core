#ifndef PPA_ADAPT_H
#define PPA_ADAPT_H

#include "pod.h"
#include <boost/fusion/include/adapt_struct.hpp>

BOOST_FUSION_ADAPT_STRUCT(
                ZLSI,
                (double, f1)
                (double, f2)
                (double, d)
                (double, sim)
                )

BOOST_FUSION_ADAPT_STRUCT(
                ZLSC,
                (double, f1)
                (double, f2)
                (FilePath, comp_solution)
                (FilePath, output_file)
                )

BOOST_FUSION_ADAPT_STRUCT(
                ZLSM,
                (FilePath, mask_file)
                (double, f2)
                (double, d)
                (double, sim)
                )

BOOST_FUSION_ADAPT_STRUCT(
                ZLSB,
                (FilePath, mask_file)
                (double, f1)
                (double, f2)
                (double, d)
                (double, sim)
                )

#endif // PPA_ADAPT_H
