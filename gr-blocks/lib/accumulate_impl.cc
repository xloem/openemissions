/* -*- c++ -*- */
/*
 * Copyright 2021 Free Software Foundation, Inc..
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <openemissions/accumulate.h>

#include <gnuradio/io_signature.h>

namespace gr {
namespace openemissions {

template <typename T>
class accumulate_impl : public accumulate<T>
{
private:
    bool d_average;
    pmt::pmt_t d_reset_tag_key;
    std::vector<T> d_accumulator;
    std::vector<T> d_counter;


public:
    /*
     * The private constructor
     */
    accumulate_impl(bool average,
                    const std::string &len_tag_key,
                    const std::string &reset_tag_key)
    : gr::tagged_stream_block("accumulate",
                              gr::io_signature::make(1, 2, sizeof(T)),
                              gr::io_signature::make(1, 2, sizeof(T)),
                              len_tag_key),
      d_average(average),
      d_reset_tag_key(pmt::string_to_symbol(reset_tag_key))
    { }
    
    /*
     * Our virtual destructor.
     */
    ~accumulate_impl()
    { }

    // Where all the action really happens
    int work(int noutput_items,
             gr_vector_int &ninput_items,
             gr_vector_const_void_star &input_items,
             gr_vector_void_star &output_items)
    {
        const T *in = reinterpret_cast<const T *>(input_items[0]);
        const T *in2 = input_items.size() > 1 ? reinterpret_cast<const T *>(input_items[1]) : nullptr;
        size_t nitems_min = in2 ? std::min(ninput_items[0], ninput_items[1]) : ninput_items[0];
        size_t nitems_max = in2 ? std::max(ninput_items[0], ninput_items[1]) : ninput_items[0];
        T *out = reinterpret_cast<T *>(output_items[0]);
        T *out2 = output_items.size() > 1 ? reinterpret_cast<T *>(output_items[1]) : nullptr;
    
        // size up accumulator if needed
        // every run we output the whole accumulator, regardless of packet size
        if (nitems_max > d_accumulator.size()) {
            d_accumulator.resize(nitems_max, 0);
            d_counter.resize(nitems_max, 0);
        }

        // grab the tags, and prepare iterators to walk them in parallel
        std::vector<std::vector<gr::tag_t>> tags(input_items.size());
        std::vector<std::vector<gr::tag_t>::iterator> tag_its(input_items.size());
        std::vector<int64_t> input_offsets(input_items.size());
        int64_t output_offset = this->nitems_written(0);
        for (size_t input = 0; input < input_items.size(); input ++) {
            this->get_tags_in_window(tags[input], input, 0, ninput_items[input]);
            tag_its[input] = tags[input].begin();
            input_offsets[input] = this->nitems_read(input);
        }

        // loop through each item of the accumulator, each packet sample
        for (size_t item = 0; item < d_accumulator.size(); item ++) {
            // if this packet has this sample, add it in
            if (item < nitems_min) {
                if (in2) {
                    d_counter[item] += in[item] * std::real(in2[item]);
                    d_accumulator[item] += std::real(in2[item]);
                } else {
                    d_accumulator[item] += in[item];
                    d_counter[item] += T(1);
                }
            }

            // output the accumulator
            out[item] = d_accumulator[item];
            if (d_average && std::real(d_counter[item])) {
                out[item] /= std::real(d_counter[item]);
            }

            // output the counter
            if (out2) {
                out2[item] = d_counter[item];
            }

            // process the tags to match where we are in iterating the accumulator
            for (size_t input = 0; input < input_items.size(); input ++) {
                uint64_t item_offset = input_offsets[input] + item;
                for (
                    std::vector<gr::tag_t>::iterator & tag_it = tag_its[input];
                    tag_it != tags[input].end() && tag_it->offset <= item_offset;
                    tag_it ++
                ) {
                    gr::tag_t tag = *tag_it;
                    // if found, the reset tag empties or scales the accumulator
                    if (tag.key == d_reset_tag_key) {
                        double value = 0;
                        if (!pmt::is_null(tag.value)) {
                            value = pmt::to_double(tag.value);
                        }
                        for (size_t j = 0; j < d_counter.size(); j ++) {
                            d_accumulator[j] *= value;
                            d_counter[j] *= value;
                        }
                    }
                }
            }
        }

        return d_accumulator.size();
    }


protected:
    int calculate_output_stream_length(const gr_vector_int &ninput_items) {
        size_t length = std::max((size_t)ninput_items[0], d_accumulator.size());
        if (ninput_items.size() > 1) {
            length = std::max(length, (size_t)ninput_items[1]);
        }
        return length;
    }
};

template <typename T>
typename accumulate<T>::sptr
accumulate<T>::make(bool average,
                    const std::string &len_tag_key,
                    const std::string &reset_tag_key)
{
    return gnuradio::make_block_sptr<accumulate_impl<T>>(average, len_tag_key, reset_tag_key);
}

template class accumulate<std::complex<double>>;
template class accumulate<double>;
template class accumulate<std::complex<float>>;
template class accumulate<float>;
template class accumulate<uint64_t>;
template class accumulate<int64_t>;
template class accumulate<int>;

} /* namespace openemissions */
} /* namespace gr */
