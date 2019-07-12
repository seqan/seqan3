// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/contrib/stream/gz_istream.hpp>
#include <seqan3/io/alignment_file/format_bam.hpp>
#include <seqan3/range/decorator/gap_decorator.hpp>

#include "alignment_file_format_test_template.hpp"

template <>
struct alignment_file_read<format_bam> : public alignment_file_data
{
    // -----------------------------------------------------------------------------------------------------------------
    // formatted input
    // -----------------------------------------------------------------------------------------------------------------
    // See format_sam_test for the corresponding input in human readable-form.
    // The byte sequence here are all uncompressed bam files since the file handles compression.

    using stream_type = std::istringstream;

    std::string big_header_input{
        '\x42', '\x41', '\x4D', '\x01', '\xB7', '\x01', '\x00', '\x00', '\x40', '\x48', '\x44', '\x09', '\x56',
        '\x4E', '\x3A', '\x31', '\x2E', '\x36', '\x09', '\x53', '\x4F', '\x3A', '\x63', '\x6F', '\x6F', '\x72',
        '\x64', '\x69', '\x6E', '\x61', '\x74', '\x65', '\x09', '\x53', '\x53', '\x3A', '\x63', '\x6F', '\x6F',
        '\x72', '\x64', '\x69', '\x6E', '\x61', '\x74', '\x65', '\x3A', '\x71', '\x75', '\x65', '\x72', '\x79',
        '\x6E', '\x61', '\x6D', '\x65', '\x09', '\x47', '\x4F', '\x3A', '\x6E', '\x6F', '\x6E', '\x65', '\x0A',
        '\x40', '\x50', '\x47', '\x09', '\x49', '\x44', '\x3A', '\x71', '\x63', '\x09', '\x50', '\x4E', '\x3A',
        '\x71', '\x75', '\x61', '\x6C', '\x69', '\x74', '\x79', '\x5F', '\x63', '\x6F', '\x6E', '\x74', '\x72',
        '\x6F', '\x6C', '\x09', '\x43', '\x4C', '\x3A', '\x71', '\x63', '\x20', '\x2D', '\x66', '\x20', '\x66',
        '\x69', '\x6C', '\x65', '\x31', '\x09', '\x44', '\x53', '\x3A', '\x74', '\x72', '\x69', '\x6D', '\x20',
        '\x72', '\x65', '\x61', '\x64', '\x73', '\x20', '\x77', '\x69', '\x74', '\x68', '\x20', '\x6C', '\x6F',
        '\x77', '\x20', '\x71', '\x75', '\x61', '\x6C', '\x09', '\x56', '\x4E', '\x3A', '\x31', '\x2E', '\x30',
        '\x2E', '\x30', '\x0A', '\x40', '\x50', '\x47', '\x09', '\x49', '\x44', '\x3A', '\x6E', '\x6F', '\x76',
        '\x6F', '\x61', '\x6C', '\x69', '\x67', '\x6E', '\x09', '\x50', '\x4E', '\x3A', '\x6E', '\x6F', '\x76',
        '\x6F', '\x61', '\x6C', '\x69', '\x67', '\x6E', '\x09', '\x56', '\x4E', '\x3A', '\x56', '\x33', '\x2E',
        '\x30', '\x32', '\x2E', '\x30', '\x37', '\x09', '\x43', '\x4C', '\x3A', '\x6E', '\x6F', '\x76', '\x6F',
        '\x61', '\x6C', '\x69', '\x67', '\x6E', '\x20', '\x2D', '\x64', '\x20', '\x2F', '\x70', '\x61', '\x74',
        '\x68', '\x2F', '\x68', '\x73', '\x33', '\x37', '\x64', '\x35', '\x2E', '\x6E', '\x64', '\x78', '\x20',
        '\x2D', '\x66', '\x20', '\x2F', '\x70', '\x61', '\x74', '\x68', '\x2F', '\x66', '\x69', '\x6C', '\x65',
        '\x2E', '\x66', '\x61', '\x73', '\x74', '\x71', '\x2E', '\x67', '\x7A', '\x09', '\x50', '\x50', '\x3A',
        '\x71', '\x63', '\x0A', '\x40', '\x53', '\x51', '\x09', '\x53', '\x4E', '\x3A', '\x72', '\x65', '\x66',
        '\x09', '\x4C', '\x4E', '\x3A', '\x32', '\x34', '\x39', '\x32', '\x35', '\x30', '\x36', '\x32', '\x31',
        '\x0A', '\x40', '\x53', '\x51', '\x09', '\x53', '\x4E', '\x3A', '\x72', '\x65', '\x66', '\x32', '\x09',
        '\x4C', '\x4E', '\x3A', '\x32', '\x34', '\x33', '\x31', '\x39', '\x39', '\x33', '\x37', '\x33', '\x09',
        '\x41', '\x53', '\x3A', '\x68', '\x73', '\x33', '\x37', '\x64', '\x35', '\x0A', '\x40', '\x52', '\x47',
        '\x09', '\x49', '\x44', '\x3A', '\x55', '\x30', '\x61', '\x5F', '\x41', '\x32', '\x5F', '\x4C', '\x31',
        '\x09', '\x50', '\x4C', '\x3A', '\x69', '\x6C', '\x6C', '\x75', '\x6D', '\x69', '\x6E', '\x61', '\x09',
        '\x50', '\x55', '\x3A', '\x31', '\x09', '\x4C', '\x42', '\x3A', '\x31', '\x09', '\x53', '\x4D', '\x3A',
        '\x4E', '\x41', '\x31', '\x32', '\x38', '\x37', '\x38', '\x0A', '\x40', '\x52', '\x47', '\x09', '\x49',
        '\x44', '\x3A', '\x55', '\x30', '\x61', '\x5F', '\x41', '\x32', '\x5F', '\x4C', '\x32', '\x09', '\x50',
        '\x4C', '\x3A', '\x69', '\x6C', '\x6C', '\x75', '\x6D', '\x69', '\x6E', '\x61', '\x09', '\x53', '\x4D',
        '\x3A', '\x4E', '\x41', '\x31', '\x32', '\x38', '\x37', '\x38', '\x09', '\x50', '\x55', '\x3A', '\x31',
        '\x09', '\x4C', '\x42', '\x3A', '\x31', '\x0A', '\x40', '\x43', '\x4F', '\x09', '\x54', '\x72', '\x61',
        '\x6C', '\x61', '\x6C', '\x61', '\x6C', '\x61', '\x6C', '\x61', '\x6C', '\x61', '\x6C', '\x61', '\x20',
        '\x74', '\x68', '\x69', '\x73', '\x20', '\x69', '\x73', '\x20', '\x61', '\x20', '\x63', '\x6F', '\x6D',
        '\x6D', '\x65', '\x6E', '\x74', '\x0A', '\x02', '\x00', '\x00', '\x00', '\x04', '\x00', '\x00', '\x00',
        '\x72', '\x65', '\x66', '\x00', '\x3D', '\x43', '\xDB', '\x0E', '\x05', '\x00', '\x00', '\x00', '\x72',
        '\x65', '\x66', '\x32', '\x00', '\x8D', '\xED', '\x7E', '\x0E'
    };

    std::string simple_three_reads_input{
        '\x42', '\x41', '\x4D', '\x01', '\x1C', '\x00', '\x00', '\x00', '\x40', '\x48', '\x44', '\x09', '\x56',
        '\x4E', '\x3A', '\x31', '\x2E', '\x36', '\x0A', '\x40', '\x53', '\x51', '\x09', '\x53', '\x4E', '\x3A',
        '\x72', '\x65', '\x66', '\x09', '\x4C', '\x4E', '\x3A', '\x33', '\x34', '\x0A', '\x01', '\x00', '\x00',
        '\x00', '\x04', '\x00', '\x00', '\x00', '\x72', '\x65', '\x66', '\x00', '\x22', '\x00', '\x00', '\x00',
        '\x48', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x06',
        '\x3D', '\x49', '\x12', '\x05', '\x00', '\x29', '\x00', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00',
        '\x00', '\x00', '\x09', '\x00', '\x00', '\x00', '\x2C', '\x01', '\x00', '\x00', '\x72', '\x65', '\x61',
        '\x64', '\x31', '\x00', '\x14', '\x00', '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x12', '\x00',
        '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x11', '\x00', '\x00', '\x00', '\x12', '\x48', '\x00',
        '\x02', '\x02', '\x03', '\x41', '\x53', '\x43', '\x02', '\x4E', '\x4D', '\x43', '\x07', '\x56', '\x00',
        '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x01', '\x00', '\x00', '\x00', '\x06', '\x3E', '\x49',
        '\x12', '\x05', '\x00', '\x2A', '\x00', '\x09', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00',
        '\x09', '\x00', '\x00', '\x00', '\x2C', '\x01', '\x00', '\x00', '\x72', '\x65', '\x61', '\x64', '\x32',
        '\x00', '\x15', '\x00', '\x00', '\x00', '\x70', '\x00', '\x00', '\x00', '\x12', '\x00', '\x00', '\x00',
        '\x10', '\x00', '\x00', '\x00', '\x14', '\x00', '\x00', '\x00', '\x14', '\x42', '\x84', '\xF1', '\x40',
        '\x00', '\x02', '\x02', '\x03', '\x05', '\x06', '\x07', '\x08', '\x09', '\x78', '\x79', '\x42', '\x53',
        '\x03', '\x00', '\x00', '\x00', '\x03', '\x00', '\x04', '\x00', '\x05', '\x00', '\x5A', '\x00', '\x00',
        '\x00', '\x00', '\x00', '\x00', '\x00', '\x02', '\x00', '\x00', '\x00', '\x06', '\x3F', '\x49', '\x12',
        '\x0A', '\x00', '\x2B', '\x00', '\x08', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x09',
        '\x00', '\x00', '\x00', '\x2C', '\x01', '\x00', '\x00', '\x72', '\x65', '\x61', '\x64', '\x33', '\x00',
        '\x14', '\x00', '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x12', '\x00', '\x00', '\x00', '\x10',
        '\x00', '\x00', '\x00', '\x11', '\x00', '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x11', '\x00',
        '\x00', '\x00', '\x12', '\x00', '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x14', '\x00', '\x00',
        '\x00', '\x44', '\x14', '\x81', '\x81', '\x00', '\x00', '\x09', '\x0A', '\x0B', '\x0C', '\x0D', '\x0E'
    };

    std::string verbose_reads_input{
        '\x42', '\x41', '\x4D', '\x01', '\xA3', '\x00', '\x00', '\x00', '\x40', '\x48', '\x44', '\x09', '\x56',
        '\x4E', '\x3A', '\x31', '\x2E', '\x36', '\x09', '\x53', '\x4F', '\x3A', '\x75', '\x6E', '\x6B', '\x6E',
        '\x6F', '\x77', '\x6E', '\x09', '\x47', '\x4F', '\x3A', '\x6E', '\x6F', '\x6E', '\x65', '\x0A', '\x40',
        '\x53', '\x51', '\x09', '\x53', '\x4E', '\x3A', '\x72', '\x65', '\x66', '\x09', '\x4C', '\x4E', '\x3A',
        '\x33', '\x34', '\x09', '\x41', '\x4E', '\x3A', '\x6F', '\x74', '\x68', '\x65', '\x72', '\x5F', '\x6E',
        '\x61', '\x6D', '\x65', '\x0A', '\x40', '\x52', '\x47', '\x09', '\x49', '\x44', '\x3A', '\x67', '\x72',
        '\x6F', '\x75', '\x70', '\x31', '\x09', '\x6D', '\x6F', '\x72', '\x65', '\x20', '\x69', '\x6E', '\x66',
        '\x6F', '\x0A', '\x40', '\x50', '\x47', '\x09', '\x49', '\x44', '\x3A', '\x70', '\x72', '\x6F', '\x67',
        '\x31', '\x09', '\x50', '\x4E', '\x3A', '\x63', '\x6F', '\x6F', '\x6C', '\x5F', '\x70', '\x72', '\x6F',
        '\x67', '\x72', '\x61', '\x6D', '\x09', '\x43', '\x4C', '\x3A', '\x2E', '\x2F', '\x70', '\x72', '\x6F',
        '\x67', '\x31', '\x09', '\x50', '\x50', '\x3A', '\x61', '\x09', '\x44', '\x53', '\x3A', '\x62', '\x09',
        '\x56', '\x4E', '\x3A', '\x63', '\x0A', '\x40', '\x43', '\x4F', '\x09', '\x54', '\x68', '\x69', '\x73',
        '\x20', '\x69', '\x73', '\x20', '\x61', '\x20', '\x63', '\x6F', '\x6D', '\x6D', '\x65', '\x6E', '\x74',
        '\x2E', '\x0A', '\x01', '\x00', '\x00', '\x00', '\x04', '\x00', '\x00', '\x00', '\x72', '\x65', '\x66',
        '\x00', '\x22', '\x00', '\x00', '\x00', '\x64', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00',
        '\x00', '\x00', '\x00', '\x00', '\x06', '\x3D', '\x49', '\x12', '\x05', '\x00', '\x29', '\x00', '\x04',
        '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x09', '\x00', '\x00', '\x00', '\x2C', '\x01',
        '\x00', '\x00', '\x72', '\x65', '\x61', '\x64', '\x31', '\x00', '\x14', '\x00', '\x00', '\x00', '\x10',
        '\x00', '\x00', '\x00', '\x12', '\x00', '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x11', '\x00',
        '\x00', '\x00', '\x12', '\x48', '\x00', '\x02', '\x02', '\x03', '\x41', '\x53', '\x43', '\x02', '\x43',
        '\x43', '\x53', '\x2C', '\x01', '\x4E', '\x4D', '\x63', '\xF9', '\x61', '\x61', '\x41', '\x63', '\x63',
        '\x63', '\x73', '\xD4', '\xFE', '\x66', '\x66', '\x66', '\x66', '\x66', '\x46', '\x40', '\x7A', '\x7A',
        '\x5A', '\x73', '\x74', '\x72', '\x00', '\xA7', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00',
        '\x01', '\x00', '\x00', '\x00', '\x06', '\x3E', '\x49', '\x12', '\x04', '\x00', '\x2A', '\x00', '\x09',
        '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x09', '\x00', '\x00', '\x00', '\x2C', '\x01',
        '\x00', '\x00', '\x72', '\x65', '\x61', '\x64', '\x32', '\x00', '\x70', '\x00', '\x00', '\x00', '\x12',
        '\x00', '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x14', '\x00', '\x00', '\x00', '\x14', '\x42',
        '\x84', '\xF1', '\x40', '\x00', '\x02', '\x02', '\x03', '\x05', '\x06', '\x07', '\x08', '\x09', '\x62',
        '\x43', '\x42', '\x43', '\x02', '\x00', '\x00', '\x00', '\x03', '\xC8', '\x62', '\x49', '\x42', '\x49',
        '\x01', '\x00', '\x00', '\x00', '\x00', '\xD8', '\x94', '\x11', '\x62', '\x53', '\x42', '\x53', '\x03',
        '\x00', '\x00', '\x00', '\x2C', '\x01', '\x28', '\x00', '\xF4', '\x01', '\x62', '\x63', '\x42', '\x63',
        '\x01', '\x00', '\x00', '\x00', '\xFD', '\x62', '\x66', '\x42', '\x66', '\x03', '\x00', '\x00', '\x00',
        '\x00', '\x00', '\x60', '\x40', '\xCD', '\xCC', '\xCC', '\x3D', '\x33', '\x33', '\x2F', '\x42', '\x62',
        '\x69', '\x42', '\x69', '\x03', '\x00', '\x00', '\x00', '\xFD', '\xFF', '\xFF', '\xFF', '\xC8', '\x00',
        '\x00', '\x00', '\x30', '\xFE', '\xFE', '\xFF', '\x62', '\x73', '\x42', '\x73', '\x03', '\x00', '\x00',
        '\x00', '\xFD', '\xFF', '\xC8', '\x00', '\xD4', '\xFE', '\x5A', '\x00', '\x00', '\x00', '\x00', '\x00',
        '\x00', '\x00', '\x02', '\x00', '\x00', '\x00', '\x06', '\x3F', '\x49', '\x12', '\x0A', '\x00', '\x2B',
        '\x00', '\x08', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x09', '\x00', '\x00', '\x00',
        '\x2C', '\x01', '\x00', '\x00', '\x72', '\x65', '\x61', '\x64', '\x33', '\x00', '\x14', '\x00', '\x00',
        '\x00', '\x10', '\x00', '\x00', '\x00', '\x12', '\x00', '\x00', '\x00', '\x10', '\x00', '\x00', '\x00',
        '\x11', '\x00', '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x11', '\x00', '\x00', '\x00', '\x12',
        '\x00', '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x14', '\x00', '\x00', '\x00', '\x44', '\x14',
        '\x81', '\x81', '\x00', '\x00', '\x09', '\x0A', '\x0B', '\x0C', '\x0D', '\x0E'
    };

    std::string empty_input{
        '\x42', '\x41', '\x4D', '\x01', '\x1C', '\x00', '\x00', '\x00', '\x40', '\x48', '\x44', '\x09', '\x56',
        '\x4E', '\x3A', '\x31', '\x2E', '\x36', '\x0A', '\x40', '\x53', '\x51', '\x09', '\x53', '\x4E', '\x3A',
        '\x72', '\x65', '\x66', '\x09', '\x4C', '\x4E', '\x3A', '\x33', '\x34', '\x0A', '\x01', '\x00', '\x00',
        '\x00', '\x04', '\x00', '\x00', '\x00', '\x72', '\x65', '\x66', '\x00', '\x22', '\x00', '\x00', '\x00',
        '\x22', '\x00', '\x00', '\x00', '\xFF', '\xFF', '\xFF', '\xFF', '\xFF', '\xFF', '\xFF', '\xFF', '\x02',
        '\x00', '\x48', '\x12', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\xFF', '\xFF',
        '\xFF', '\xFF', '\xFF', '\xFF', '\xFF', '\xFF', '\x00', '\x00', '\x00', '\x00', '\x2A', '\x00'
    };

    std::string empty_cigar{
        '\x42', '\x41', '\x4D', '\x01', '\x1C', '\x00', '\x00', '\x00', '\x40', '\x48', '\x44', '\x09', '\x56',
        '\x4E', '\x3A', '\x31', '\x2E', '\x36', '\x0A', '\x40', '\x53', '\x51', '\x09', '\x53', '\x4E', '\x3A',
        '\x72', '\x65', '\x66', '\x09', '\x4C', '\x4E', '\x3A', '\x33', '\x34', '\x0A', '\x01', '\x00', '\x00',
        '\x00', '\x04', '\x00', '\x00', '\x00', '\x72', '\x65', '\x66', '\x00', '\x22', '\x00', '\x00', '\x00',
        '\x34', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x06',
        '\x3D', '\x49', '\x12', '\x00', '\x00', '\x2D', '\x00', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00',
        '\x00', '\x00', '\x09', '\x00', '\x00', '\x00', '\x2C', '\x01', '\x00', '\x00', '\x72', '\x65', '\x61',
        '\x64', '\x31', '\x00', '\x12', '\x48', '\x00', '\x02', '\x02', '\x03', '\x41', '\x53', '\x43', '\x02',
        '\x4E', '\x4D', '\x43', '\x07'
    };

    std::string unknown_ref{
        '\x42', '\x41', '\x4D', '\x01', '\x1C', '\x00', '\x00', '\x00', '\x40', '\x48', '\x44', '\x09', '\x56',
        '\x4E', '\x3A', '\x31', '\x2E', '\x36', '\x0A', '\x40', '\x53', '\x51', '\x09', '\x53', '\x4E', '\x3A',
        '\x72', '\x61', '\x66', '\x09', '\x4C', '\x4E', '\x3A', '\x33', '\x34', '\x0A', '\x01', '\x00', '\x00',
        '\x00', '\x04', '\x00', '\x00', '\x00', '\x72', '\x61', '\x66', '\x00', '\x22', '\x00', '\x00', '\x00',
        '\x56', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x06',
        '\x3D', '\x49', '\x12', '\x05', '\x00', '\x29', '\x00', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00',
        '\x00', '\x00', '\x09', '\x00', '\x00', '\x00', '\x2C', '\x01', '\x00', '\x00', '\x72', '\x65', '\x61',
        '\x64', '\x31', '\x00', '\x14', '\x00', '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x12', '\x00',
        '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x11', '\x00', '\x00', '\x00', '\x12', '\x48', '\x00',
        '\x02', '\x02', '\x03', '\x61', '\x61', '\x41', '\x63', '\x41', '\x53', '\x43', '\x02', '\x66', '\x66',
        '\x66', '\x66', '\x66', '\x46', '\x40', '\x7A', '\x7A', '\x5A', '\x73', '\x74', '\x72', '\x00'
    };

    /* bytes were modified to a ref id of 8448: \x00 \x00 \x21 \x00*/
    std::string unknown_ref_header{
        '\x42', '\x41', '\x4D', '\x01', '\x1C', '\x00', '\x00', '\x00', '\x40', '\x48', '\x44', '\x09', '\x56',
        '\x4E', '\x3A', '\x31', '\x2E', '\x36', '\x0A', '\x40', '\x53', '\x51', '\x09', '\x53', '\x4E', '\x3A',
        '\x72', '\x65', '\x66', '\x09', '\x4C', '\x4E', '\x3A', '\x33', '\x34', '\x0A', '\x01', '\x00', '\x00',
        '\x00', '\x04', '\x00', '\x00', '\x00', '\x72', '\x65', '\x66', '\x00', '\x22', '\x00', '\x00', '\x00',
        '\x56', '\x00', '\x00', '\x00', '\x00', '\x21', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x06',
        '\x3D', '\x49', '\x12', '\x05', '\x00', '\x29', '\x00', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00',
        '\x00', '\x00', '\x09', '\x00', '\x00', '\x00', '\x2C', '\x01', '\x00', '\x00', '\x72', '\x65', '\x61',
        '\x64', '\x31', '\x00', '\x14', '\x00', '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x12', '\x00',
        '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x11', '\x00', '\x00', '\x00', '\x12', '\x48', '\x00',
        '\x02', '\x02', '\x03', '\x61', '\x61', '\x41', '\x63', '\x41', '\x53', '\x43', '\x02', '\x66', '\x66',
        '\x66', '\x66', '\x66', '\x46', '\x40', '\x7A', '\x7A', '\x5A', '\x73', '\x74', '\x72', '\x00', '\x0A'
    };

    std::string simple_three_reads_output{ // no hard clipping in output
        '\x42', '\x41', '\x4D', '\x01', '\x1C', '\x00', '\x00', '\x00', '\x40', '\x48', '\x44', '\x09', '\x56',
        '\x4E', '\x3A', '\x31', '\x2E', '\x36', '\x0A', '\x40', '\x53', '\x51', '\x09', '\x53', '\x4E', '\x3A',
        '\x72', '\x65', '\x66', '\x09', '\x4C', '\x4E', '\x3A', '\x33', '\x34', '\x0A', '\x01', '\x00', '\x00',
        '\x00', '\x04', '\x00', '\x00', '\x00', '\x72', '\x65', '\x66', '\x00', '\x22', '\x00', '\x00', '\x00',
        '\x48', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x06',
        '\x3D', '\x49', '\x12', '\x05', '\x00', '\x29', '\x00', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00',
        '\x00', '\x00', '\x09', '\x00', '\x00', '\x00', '\x2C', '\x01', '\x00', '\x00', '\x72', '\x65', '\x61',
        '\x64', '\x31', '\x00', '\x14', '\x00', '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x12', '\x00',
        '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x11', '\x00', '\x00', '\x00', '\x12', '\x48', '\x00',
        '\x02', '\x02', '\x03', '\x41', '\x53', '\x43', '\x02', '\x4E', '\x4D', '\x43', '\x07', '\x52', '\x00',
        '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x01', '\x00', '\x00', '\x00', '\x06', '\x3E', '\x49',
        '\x12', '\x04', '\x00', '\x2A', '\x00', '\x09', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00',
        '\x09', '\x00', '\x00', '\x00', '\x2C', '\x01', '\x00', '\x00', '\x72', '\x65', '\x61', '\x64', '\x32',
        '\x00', '\x70', '\x00', '\x00', '\x00', '\x12', '\x00', '\x00', '\x00', '\x10', '\x00', '\x00', '\x00',
        '\x14', '\x00', '\x00', '\x00', '\x14', '\x42', '\x84', '\xF1', '\x40', '\x00', '\x02', '\x02', '\x03',
        '\x05', '\x06', '\x07', '\x08', '\x09', '\x78', '\x79', '\x42', '\x53', '\x03', '\x00', '\x00', '\x00',
        '\x03', '\x00', '\x04', '\x00', '\x05', '\x00', '\x5A', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00',
        '\x00', '\x02', '\x00', '\x00', '\x00', '\x06', '\x3F', '\x49', '\x12', '\x0A', '\x00', '\x2B', '\x00',
        '\x08', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x09', '\x00', '\x00', '\x00', '\x2C',
        '\x01', '\x00', '\x00', '\x72', '\x65', '\x61', '\x64', '\x33', '\x00', '\x14', '\x00', '\x00', '\x00',
        '\x10', '\x00', '\x00', '\x00', '\x12', '\x00', '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x11',
        '\x00', '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x11', '\x00', '\x00', '\x00', '\x12', '\x00',
        '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x14', '\x00', '\x00', '\x00', '\x44', '\x14', '\x81',
        '\x81', '\x00', '\x00', '\x09', '\x0A', '\x0B', '\x0C', '\x0D', '\x0E'
    };

    std::string verbose_output{verbose_reads_input};

    std::string special_output{
        '\x42', '\x41', '\x4D', '\x01', '\x1C', '\x00', '\x00', '\x00', '\x40', '\x48', '\x44', '\x09', '\x56',
        '\x4E', '\x3A', '\x31', '\x2E', '\x36', '\x0A', '\x40', '\x53', '\x51', '\x09', '\x53', '\x4E', '\x3A',
        '\x72', '\x65', '\x66', '\x09', '\x4C', '\x4E', '\x3A', '\x33', '\x34', '\x0A', '\x01', '\x00', '\x00',
        '\x00', '\x04', '\x00', '\x00', '\x00', '\x72', '\x65', '\x66', '\x00', '\x22', '\x00', '\x00', '\x00',
        '\x40', '\x00', '\x00', '\x00', '\xFF', '\xFF', '\xFF', '\xFF', '\x00', '\x00', '\x00', '\x00', '\x06',
        '\x3D', '\x49', '\x12', '\x05', '\x00', '\x29', '\x00', '\x04', '\x00', '\x00', '\x00', '\xFF', '\xFF',
        '\xFF', '\xFF', '\xFF', '\xFF', '\xFF', '\xFF', '\x00', '\x00', '\x00', '\x00', '\x72', '\x65', '\x61',
        '\x64', '\x31', '\x00', '\x14', '\x00', '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x12', '\x00',
        '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x11', '\x00', '\x00', '\x00', '\x12', '\x48', '\x00',
        '\x02', '\x02', '\x03'
    };
};

// ---------------------------------------------------------------------------------------------------------------------
// parametrized tests
// ---------------------------------------------------------------------------------------------------------------------

INSTANTIATE_TYPED_TEST_CASE_P(bam, alignment_file_read, format_bam);
INSTANTIATE_TYPED_TEST_CASE_P(bam, alignment_file_write, format_bam);

// ---------------------------------------------------------------------------------------------------------------------
// BAM specifics
// ---------------------------------------------------------------------------------------------------------------------

struct bam_format : public alignment_file_data
{};

TEST_F(bam_format, wrong_magic_bytes)
{
    std::string cam1_bute_str{'\x43', '\x41', '\x4D', '\x01' /*CAM\1*/};

    std::istringstream stream{cam1_bute_str};
    detail::alignment_file_input_format<format_bam> format{};

    EXPECT_THROW(format.read(stream, input_options, std::ignore, this->header, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore),
                 format_error);

}

TEST_F(bam_format, unknown_ref_in_header)
{
    std::string unknown_ref{ // raf instead of ref
        // @HD     VN:1.0
        // @SQ     SN:raf  LN:34
        '\x42', '\x41', '\x4D', '\x01', '\x1C', '\x00', '\x00', '\x00', '\x40', '\x48', '\x44', '\x09', '\x56',
        '\x4E', '\x3A', '\x31', '\x2E', '\x36', '\x0A', '\x40', '\x53', '\x51', '\x09', '\x53', '\x4E', '\x3A',
        '\x72', '\x65', '\x66', '\x09', '\x4C', '\x4E', '\x3A', '\x33', '\x34', '\x0A', '\x01', '\x00', '\x00',
        '\x00', '\x04', '\x00', '\x00', '\x00', '\x72', '\x61', '\x66', '\x00', '\x22', '\x00', '\x00', '\x00',
    };

    std::istringstream stream{unknown_ref};
    detail::alignment_file_input_format<format_bam> format{};

    EXPECT_THROW(format.read(stream, input_options, this->ref_sequences, this->header, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore),
                 format_error);
}

TEST_F(bam_format, wrong_ref_length_in_header)
{
    std::string wrong_ref_length{ // 33 instead of 34
        // @HD     VN:1.0
        // @SQ     SN:ref  LN:33
        '\x42', '\x41', '\x4D', '\x01', '\x1C', '\x00', '\x00', '\x00', '\x40', '\x48', '\x44', '\x09', '\x56',
        '\x4E', '\x3A', '\x31', '\x2E', '\x36', '\x0A', '\x40', '\x53', '\x51', '\x09', '\x53', '\x4E', '\x3A',
        '\x72', '\x65', '\x66', '\x09', '\x4C', '\x4E', '\x3A', '\x33', '\x34', '\x0A', '\x01', '\x00', '\x00',
        '\x00', '\x04', '\x00', '\x00', '\x00', '\x72', '\x65', '\x66', '\x00', '\x23', '\x00', '\x00', '\x00',
    };

    std::istringstream stream{wrong_ref_length};
    detail::alignment_file_input_format<format_bam> format{};

    EXPECT_THROW(format.read(stream, input_options, this->ref_sequences, this->header, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore),
                 format_error);
}

TEST_F(bam_format, wrong_order_in_header)
{
    std::vector<std::string> rids = {"ref", "raf"};
    alignment_file_header hdr{rids};
    hdr.ref_id_info.emplace_back(34, "");
    hdr.ref_id_info.emplace_back(30, "");
    hdr.ref_dict[hdr.ref_ids()[0]] = 0;
    hdr.ref_dict[hdr.ref_ids()[1]] = 1;

    std::string wrong_order{ // raf is first in file but second in hdr
        // @HD     VN:1.6
        // @SQ     SN:raf  LN:30
        // @SQ     SN:ref  LN:34
        '\x42', '\x41', '\x4D', '\x01', '\x2D', '\x00', '\x00', '\x00', '\x40', '\x48', '\x44', '\x09', '\x56',
        '\x4E', '\x3A', '\x31', '\x2E', '\x36', '\x0A', '\x40', '\x53', '\x51', '\x09', '\x53', '\x4E', '\x3A',
        '\x72', '\x61', '\x66', '\x09', '\x4C', '\x4E', '\x3A', '\x33', '\x30', '\x0A', '\x40', '\x53', '\x51',
        '\x09', '\x53', '\x4E', '\x3A', '\x72', '\x65', '\x66', '\x09', '\x4C', '\x4E', '\x3A', '\x33', '\x34',
        '\x0A', '\x02', '\x00', '\x00', '\x00', '\x04', '\x00', '\x00', '\x00', '\x72', '\x61', '\x66', '\x00',
        '\x1E', '\x00', '\x00', '\x00', '\x04', '\x00', '\x00', '\x00', '\x72', '\x65', '\x66', '\x00', '\x22',
        '\x00', '\x00', '\x00'
    };

    std::istringstream stream{wrong_order};
    detail::alignment_file_input_format<format_bam> format{};

    EXPECT_THROW(format.read(stream, input_options, this->ref_sequences, hdr, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore),
                 format_error);
}

TEST_F(bam_format, wrong_char_as_tag_identifier)
{
    dna5_vector seq;
    std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>> alignment;
    sam_tag_dictionary tag_dict;

    {
        std::string wrong_char_in_tag{ // Y in CG tag
            // @HD     VN:1.0
            // @SQ     SN:ref  LN:34
            // read1   41      ref     1       61      4S3N    =       10      300     ACGT    !##$    CG:Y:1S1M1D1M1I
            '\x42', '\x41', '\x4D', '\x01', '\x1C', '\x00', '\x00', '\x00', '\x40', '\x48', '\x44', '\x09', '\x56',
            '\x4E', '\x3A', '\x31', '\x2E', '\x36', '\x0A', '\x40', '\x53', '\x51', '\x09', '\x53', '\x4E', '\x3A',
            '\x72', '\x65', '\x66', '\x09', '\x4C', '\x4E', '\x3A', '\x33', '\x34', '\x0A', '\x01', '\x00', '\x00',
            '\x00', '\x04', '\x00', '\x00', '\x00', '\x72', '\x65', '\x66', '\x00', '\x22', '\x00', '\x00', '\x00',
            '\x42', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x06',
            '\x3D', '\x49', '\x12', '\x02', '\x00', '\x29', '\x00', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00',
            '\x00', '\x00', '\x09', '\x00', '\x00', '\x00', '\x2C', '\x01', '\x00', '\x00', '\x72', '\x65', '\x61',
            '\x64', '\x31', '\x00', '\x44', '\x00', '\x00', '\x00', '\x33', '\x00', '\x00', '\x00', '\x12', '\x48',
            '\x00', '\x02', '\x02', '\x03', '\x43', '\x47', '\x59', '\x31', '\x53', '\x31', '\x4D', '\x31', '\x44',
            '\x31', '\x4D', '\x31', '\x49', '\x00'
        };

        std::istringstream stream{wrong_char_in_tag};
        detail::alignment_file_input_format<format_bam> format{};

        EXPECT_THROW(format.read(stream, input_options, this->ref_sequences, this->header, seq, std::ignore,
                                    std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, alignment, std::ignore,
                                    std::ignore, std::ignore, tag_dict, std::ignore, std::ignore),
                    format_error);
    }
    {
        std::string wrong_char_in_tag{ // Y in CG:B array tag
            // @HD     VN:1.0
            // @SQ     SN:ref  LN:34
            // read1   41      ref     1       61      4S3N    =       10      300     ACGT    !##$    CG:B:Y1S1M1D1M1
            '\x42', '\x41', '\x4D', '\x01', '\x1C', '\x00', '\x00', '\x00', '\x40', '\x48', '\x44', '\x09', '\x56',
            '\x4E', '\x3A', '\x31', '\x2E', '\x36', '\x0A', '\x40', '\x53', '\x51', '\x09', '\x53', '\x4E', '\x3A',
            '\x72', '\x65', '\x66', '\x09', '\x4C', '\x4E', '\x3A', '\x33', '\x34', '\x0A', '\x01', '\x00', '\x00',
            '\x00', '\x04', '\x00', '\x00', '\x00', '\x72', '\x65', '\x66', '\x00', '\x22', '\x00', '\x00', '\x00',
            '\x42', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x06',
            '\x3D', '\x49', '\x12', '\x02', '\x00', '\x29', '\x00', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00',
            '\x00', '\x00', '\x09', '\x00', '\x00', '\x00', '\x2C', '\x01', '\x00', '\x00', '\x72', '\x65', '\x61',
            '\x64', '\x31', '\x00', '\x44', '\x00', '\x00', '\x00', '\x33', '\x00', '\x00', '\x00', '\x12', '\x48',
            '\x00', '\x02', '\x02', '\x03', '\x43', '\x47', '\x42', '\x59', '\x53', '\x31', '\x4D', '\x31', '\x44',
            '\x31', '\x4D', '\x31', '\x49', '\x00'
        };

        std::istringstream stream{wrong_char_in_tag};
        detail::alignment_file_input_format<format_bam> format{};

        EXPECT_THROW(format.read(stream, input_options, this->ref_sequences, this->header, seq, std::ignore,
                                    std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, alignment, std::ignore,
                                    std::ignore, std::ignore, tag_dict, std::ignore, std::ignore),
                    format_error);
    }
}

TEST_F(bam_format, too_long_cigar_string_read)
{
    std::string sam_file_with_too_long_cigar_string{
        // @HD     VN:1.0
        // @SQ     SN:ref  LN:34
        // read1   41      ref     1       61      4S3N    =       10      300     ACGT    !##$    CG:Z:1S1M1D1M1I
        '\x42', '\x41', '\x4D', '\x01', '\x1C', '\x00', '\x00', '\x00', '\x40', '\x48', '\x44', '\x09', '\x56',
        '\x4E', '\x3A', '\x31', '\x2E', '\x36', '\x0A', '\x40', '\x53', '\x51', '\x09', '\x53', '\x4E', '\x3A',
        '\x72', '\x65', '\x66', '\x09', '\x4C', '\x4E', '\x3A', '\x33', '\x34', '\x0A', '\x01', '\x00', '\x00',
        '\x00', '\x04', '\x00', '\x00', '\x00', '\x72', '\x65', '\x66', '\x00', '\x22', '\x00', '\x00', '\x00',
        '\x42', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x06',
        '\x3D', '\x49', '\x12', '\x02', '\x00', '\x29', '\x00', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00',
        '\x00', '\x00', '\x09', '\x00', '\x00', '\x00', '\x2C', '\x01', '\x00', '\x00', '\x72', '\x65', '\x61',
        '\x64', '\x31', '\x00', '\x44', '\x00', '\x00', '\x00', '\x33', '\x00', '\x00', '\x00', '\x12', '\x48',
        '\x00', '\x02', '\x02', '\x03', '\x43', '\x47', '\x5A', '\x31', '\x53', '\x31', '\x4D', '\x31', '\x44',
        '\x31', '\x4D', '\x31', '\x49', '\x00'
    };

    dna5_vector seq;
    std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>> alignment;
    sam_tag_dictionary tag_dict;

    {   // successful reading
        std::istringstream stream{sam_file_with_too_long_cigar_string};

        detail::alignment_file_input_format<format_bam> format{};

        ASSERT_NO_THROW(format.read(stream, input_options, this->ref_sequences, this->header, seq, std::ignore,
                                    std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, alignment, std::ignore,
                                    std::ignore, std::ignore, tag_dict, std::ignore, std::ignore));

        EXPECT_EQ(get<0>(alignment), get<0>(this->alignments[0]));
        EXPECT_EQ(get<1>(alignment), get<1>(this->alignments[0]));
        EXPECT_EQ(tag_dict.size(), 0u); // redundant CG tag is removed
    }

    {   // error: sam_tag_dictionary is not read
        std::istringstream stream{sam_file_with_too_long_cigar_string};

        detail::alignment_file_input_format<format_bam> format{};

        ASSERT_THROW(format.read(stream, input_options, this->ref_sequences, this->header, seq, std::ignore,
                                 std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, alignment, std::ignore,
                                 std::ignore, std::ignore, std::ignore, std::ignore, std::ignore), format_error);
    }

    {   // error: sequence is not read
        std::istringstream stream{sam_file_with_too_long_cigar_string};

        detail::alignment_file_input_format<format_bam> format{};

        ASSERT_THROW(format.read(stream, input_options, this->ref_sequences, this->header, std::ignore, std::ignore,
                                 std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, alignment, std::ignore,
                                 std::ignore, std::ignore, tag_dict, std::ignore, std::ignore), format_error);
    }

    {   // error no CG tag
        std::istringstream stream{std::string{
            // @HD     VN:1.0
            // @SQ     SN:ref  LN:34
            // read1   41      ref     1       61      4S3N    =       10      300     ACGT    !##$
            '\x42', '\x41', '\x4D', '\x01', '\x1C', '\x00', '\x00', '\x00', '\x40', '\x48', '\x44', '\x09', '\x56',
            '\x4E', '\x3A', '\x31', '\x2E', '\x36', '\x0A', '\x40', '\x53', '\x51', '\x09', '\x53', '\x4E', '\x3A',
            '\x72', '\x65', '\x66', '\x09', '\x4C', '\x4E', '\x3A', '\x33', '\x34', '\x0A', '\x01', '\x00', '\x00',
            '\x00', '\x04', '\x00', '\x00', '\x00', '\x72', '\x65', '\x66', '\x00', '\x22', '\x00', '\x00', '\x00',
            '\x34', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x06',
            '\x3D', '\x49', '\x12', '\x02', '\x00', '\x29', '\x00', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00',
            '\x00', '\x00', '\x09', '\x00', '\x00', '\x00', '\x2C', '\x01', '\x00', '\x00', '\x72', '\x65', '\x61',
            '\x64', '\x31', '\x00', '\x44', '\x00', '\x00', '\x00', '\x33', '\x00', '\x00', '\x00', '\x12', '\x48',
            '\x00', '\x02', '\x02', '\x03'
        }};

        detail::alignment_file_input_format<format_bam> format{};
        tag_dict.clear();

        ASSERT_THROW(format.read(stream, input_options, this->ref_sequences, this->header, seq, std::ignore,
                                 std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, alignment, std::ignore,
                                 std::ignore, std::ignore, tag_dict, std::ignore, std::ignore), format_error);
    }
}

TEST_F(bam_format, too_long_cigar_string_write)
{
    // create an alignment resulting more than 65535 cigar elements
    // -------------------------------------------------------------------------
    auto read = view::repeat_n('T'_dna5,  70'000);
    auto ref  = view::repeat_n('A'_dna5, 2 * read.size() - 1);

    auto gapped_ref  = gap_decorator{ref};

    // a gap_decorator on a repeat_n view also works but is slow when inserting gaps.
    std::vector<gapped<dna5>> gapped_read;
    gapped_read.reserve(2 * read.size());
    // create gap of length one every second character => T-T-T-T-T-T...
    for (auto chr : read)
        gapped_read.push_back(chr), gapped_read.push_back(gap{});
    gapped_read.pop_back(); // remove last gap

    auto alignment = std::tie(gapped_ref, gapped_read);

    // Expected output. ATTENTION this could not be validated by samtools as it does not support too long cigar strings
    // -------------------------------------------------------------------------
    std::string expected =
        std::string /*the beginning*/
        {
            '\x42', '\x41', '\x4D', '\x01', '\x20', '\x00', '\x00', '\x00', '\x40', '\x48', '\x44', '\x09', '\x56',
            '\x4E', '\x3A', '\x31', '\x2E', '\x36', '\x0A', '\x40', '\x53', '\x51', '\x09', '\x53', '\x4E', '\x3A',
            '\x72', '\x65', '\x66', '\x09', '\x4C', '\x4E', '\x3A', '\x31', '\x33', '\x39', '\x39', '\x39', '\x39',
            '\x0A', '\x01', '\x00', '\x00', '\x00', '\x04', '\x00', '\x00', '\x00', '\x72', '\x65', '\x66', '\x00',
            '\xDF', '\x22', '\x02', '\x00', '\x1C', '\xE0', '\x05', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00',
            '\x00', '\x00', '\x00', '\x0A', '\xFF', '\x49', '\x00', '\x02', '\x00', '\x00', '\x00', '\x70', '\x11',
            '\x01', '\x00', '\xFF', '\xFF', '\xFF', '\xFF', '\xFF', '\xFF', '\xFF', '\xFF', '\x00', '\x00', '\x00',
            '\x00', '\x6C', '\x6F', '\x6E', '\x67', '\x5F', '\x72', '\x65', '\x61', '\x64', '\x00', '\x04', '\x17',
            '\x11', '\x00', '\xF3', '\x2D', '\x22', '\x00'
        } +
        std::string(view::repeat_n('\x88', (read.size() + 1) / 2))  /*seq */ +
        std::string(view::repeat_n('\xFF',  read.size()))           /*qual*/ +
        std::string /*the beginning*/
        {
            '\x43', '\x47', '\x5A' /*tag info: CGZ*/
        };

    for (size_t i = 0; i < read.size() - 1; ++i)
        expected.append("\x31\x4D\x31\x44"/*1M1D*/);
    expected.append("\x31\x4D"/*1M*/);
    expected.push_back('\x00' /*\0*/);

    std::ostringstream os{};

    alignment_file_header header{std::vector<std::string>{this->ref_id}};
    header.ref_id_info.push_back({ref.size(), ""});
    header.ref_dict[this->ref_id] = 0;

    using default_mate_t  = std::tuple<std::string, std::optional<int32_t>, int32_t>;

    detail::alignment_file_output_format<format_bam> format;

    format.write(os, output_options, header, read, std::array<char, 0>{}/*empty qual*/,
                 std::string{"long_read"}, 0/*offset*/, std::string{}, 0/*ref_id*/, 0/*ref_offset*/, alignment,
                 0, 255, default_mate_t{}, sam_tag_dictionary{}, 0, 0);

    EXPECT_TRUE(os.str() == expected); // do not use EXPECT_EQ because if this fails the output will be huge :D
}
