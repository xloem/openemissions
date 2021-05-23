/* -*- c++ -*- */
/*
 * Copyright 2021 Free Software Foundation, Inc..
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_OPENEMISSIONS_ACCUMULATE_H
#define INCLUDED_OPENEMISSIONS_ACCUMULATE_H

#include <openemissions/api.h>
#include <gnuradio/tagged_stream_block.h>

namespace gr {
namespace openemissions {

/*!
 * \brief Integrate each packet with all the last
 * \ingroup openemissions
 *
 */
template <typename T>
class OPENEMISSIONS_API accumulate : virtual public gr::tagged_stream_block {
public:
    typedef std::shared_ptr<accumulate> sptr;

    /*!
     * Make an acccumulate block
     *
     * \param[in] average Average the packets by dividing by the total count.
     * \param[in] len_tag_key Name of the TSB's length tag key.
     * \param[in] reset_tag_key When seen, reset or scale stored sum by the tag's value.
     */
    static sptr make(bool average = true,
                     const std::string &len_tag_key = "packet_len",
                     const std::string &reset_tag_key = "reset_sum");
};

} // namespace openemissions
} // namespace gr

#endif /* INCLUDED_OPENEMISSIONS_ACCUMULATE_H */
