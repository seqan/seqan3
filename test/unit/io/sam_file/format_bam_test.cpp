// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/io/sam_file/format_bam.hpp>
#include <seqan3/range/decorator/gap_decorator.hpp>
#include <seqan3/range/views/to.hpp>

#include "sam_file_format_test_template.hpp"

template <>
struct sam_file_read<seqan3::format_bam> : public sam_file_data
{
    // -----------------------------------------------------------------------------------------------------------------
    // formatted input
    // -----------------------------------------------------------------------------------------------------------------
    // See format_sam_test for the corresponding input in human readable-form.
    // The byte sequence here are all uncompressed (both gzip and bgzf) bam files since the file handles compression.
    // Conversion can look like this: samtools view -u test.sam | bgzip -d
    // -u disables gzip compression, but the output is still bgzf compressed (decompression via bgzip)
    // pass --no-PG to samtools if you do not want PG tags which may be added automatically

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
        '\x02', '\x02', '\x03', '\x41', '\x53', '\x43', '\x02', '\x4E', '\x4D', '\x43', '\x07', '\x5A', '\x00',
        '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x01', '\x00', '\x00', '\x00', '\x06', '\x3E', '\x49',
        '\x12', '\x06', '\x00', '\x2A', '\x00', '\x09', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00',
        '\x09', '\x00', '\x00', '\x00', '\x2C', '\x01', '\x00', '\x00', '\x72', '\x65', '\x61', '\x64', '\x32',
        '\x00', '\x15', '\x00', '\x00', '\x00', '\x70', '\x00', '\x00', '\x00', '\x12', '\x00', '\x00', '\x00',
        '\x10', '\x00', '\x00', '\x00',
        '\x14', '\x00', '\x00', '\x00', '\x25', '\x00', '\x00', '\x00', '\x14', '\x42', '\x84', '\xF1', '\x40',
        '\x00', '\x02', '\x02', '\x03', '\x05', '\x06', '\x07', '\x08', '\x09', '\x78', '\x79', '\x42', '\x53',
        '\x03', '\x00', '\x00', '\x00', '\x03', '\x00', '\x04', '\x00', '\x05', '\x00', '\x5A', '\x00', '\x00',
        '\x00', '\x00', '\x00', '\x00', '\x00', '\x02', '\x00', '\x00', '\x00', '\x06', '\x3F', '\x49', '\x12',
        '\x0A', '\x00', '\x2B', '\x00', '\x08', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x09',
        '\x00', '\x00', '\x00', '\x2C', '\x01', '\x00', '\x00', '\x72', '\x65', '\x61', '\x64', '\x33', '\x00',
        '\x14', '\x00', '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x16', '\x00', '\x00', '\x00', '\x10',
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
        '\x00', '\x10', '\x00', '\x00', '\x00', '\x16', '\x00', '\x00', '\x00', '\x10', '\x00', '\x00', '\x00',
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
        '\x10', '\x00', '\x00', '\x00', '\x16', '\x00', '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x11',
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

INSTANTIATE_TYPED_TEST_SUITE_P(bam, sam_file_read, seqan3::format_bam, );
INSTANTIATE_TYPED_TEST_SUITE_P(bam, sam_file_write, seqan3::format_bam, );

// ---------------------------------------------------------------------------------------------------------------------
// BAM specifics
// ---------------------------------------------------------------------------------------------------------------------

struct bam_format : public sam_file_data
{};

TEST_F(bam_format, wrong_magic_bytes)
{
    std::istringstream stream{std::string{'\x43', '\x41', '\x4D', '\x01' /*CAM\1*/}};
    seqan3::sam_file_input fin{stream, seqan3::format_bam{}};
    EXPECT_THROW(fin.begin(), seqan3::format_error);
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
    seqan3::sam_file_input fin{stream, this->ref_ids, this->ref_sequences, seqan3::format_bam{}};
    EXPECT_THROW(fin.begin(), seqan3::format_error);
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
    seqan3::sam_file_input fin{stream, this->ref_ids, this->ref_sequences, seqan3::format_bam{}};
    EXPECT_THROW(fin.begin(), seqan3::format_error);
}

TEST_F(bam_format, wrong_order_in_header)
{
    std::vector<std::string> rids = {"ref", "raf"};
    std::vector<seqan3::dna5_vector> rseqs = {"ATCGAGATCGATCGATCGAGAGCTAGCGATCGAG"_dna5,
                                              "ATCGAGATCGATCGATCGAGAGCTAGCGAT"_dna5};

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
    seqan3::sam_file_input fin{stream, rids, rseqs, seqan3::format_bam{}};
    EXPECT_THROW(fin.begin(), seqan3::format_error);
}

TEST_F(bam_format, wrong_char_as_tag_identifier)
{
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
        seqan3::sam_file_input fin{stream, this->ref_ids, this->ref_sequences, seqan3::format_bam{}};
        EXPECT_THROW(fin.begin(), seqan3::format_error);
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
        seqan3::sam_file_input fin{stream, this->ref_ids, this->ref_sequences, seqan3::format_bam{}};
        EXPECT_THROW(fin.begin(), seqan3::format_error);
    }
}

TEST_F(bam_format, invalid_cigar_op)
{
    {
        std::string wrong_char_in_tag{// "1D" replaced by "1?" (D is encoded as 2, but 2 was replaced by 14)
            // @HD  VN:1.6
            // @SQ SN:ref  LN:34
            // read1   41  ref 1   61  1S1M1D1M1I  ref 10  300 ACGT    !##$    AS:i:2  NM:i:7
            '\x42', '\x41', '\x4D', '\x01', '\x1C', '\x00', '\x00', '\x00', '\x40', '\x48', '\x44', '\x09', '\x56',
            '\x4E', '\x3A', '\x31', '\x2E', '\x36', '\x0A', '\x40', '\x53', '\x51', '\x09', '\x53', '\x4E', '\x3A',
            '\x72', '\x65', '\x66', '\x09', '\x4C', '\x4E', '\x3A', '\x33', '\x34', '\x0A', '\x01', '\x00', '\x00',
            '\x00', '\x04', '\x00', '\x00', '\x00', '\x72', '\x65', '\x66', '\x00', '\x22', '\x00', '\x00', '\x00',
            '\x48', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x06',
            '\x3D', '\x49', '\x12', '\x05', '\x00', '\x29', '\x00', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00',
            '\x00', '\x00', '\x09', '\x00', '\x00', '\x00', '\x2C', '\x01', '\x00', '\x00', '\x72', '\x65', '\x61',
            '\x64', '\x31', '\x00', '\x14', '\x00', '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x1E', '\x00',
            '\x00', '\x00', '\x10', '\x00', '\x00', '\x00', '\x11', '\x00', '\x00', '\x00', '\x12', '\x48', '\x00',
            '\x02', '\x02', '\x03', '\x41', '\x53', '\x43', '\x02', '\x4E', '\x4D', '\x43', '\x07'
        };

        std::istringstream stream{wrong_char_in_tag};
        seqan3::sam_file_input fin{stream, this->ref_ids, this->ref_sequences, seqan3::format_bam{}};
        EXPECT_THROW(fin.begin(), seqan3::format_error);
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

    {   // successful reading
        std::istringstream stream{sam_file_with_too_long_cigar_string};

        seqan3::sam_file_input fin{stream, this->ref_ids, this->ref_sequences, seqan3::format_bam{}};

        EXPECT_RANGE_EQ(std::get<0>(seqan3::get<seqan3::field::alignment>(*fin.begin())),
                        std::get<0>(this->alignments[0]));
        EXPECT_RANGE_EQ(std::get<1>(seqan3::get<seqan3::field::alignment>(*fin.begin())),
                        std::get<1>(this->alignments[0]));
        EXPECT_EQ(seqan3::get<seqan3::field::tags>(*fin.begin()).size(), 0u); // redundant CG tag is removed

        EXPECT_RANGE_EQ(std::get<0>((*fin.begin()).alignment()),
                        std::get<0>(this->alignments[0]));
        EXPECT_RANGE_EQ(std::get<1>((*fin.begin()).alignment()),
                        std::get<1>(this->alignments[0]));
        EXPECT_EQ((*fin.begin()).tags().size(), 0u); // redundant CG tag is removed
    }

    {   // error: sam_tag_dictionary is not read
        std::istringstream stream{sam_file_with_too_long_cigar_string};

        seqan3::sam_file_input fin{stream, seqan3::format_bam{}, seqan3::fields<seqan3::field::alignment>{}};
        ASSERT_THROW(fin.begin(), seqan3::format_error);
    }

    {   // error: sequence is not read
        std::istringstream stream{sam_file_with_too_long_cigar_string};

        seqan3::sam_file_input fin{stream, seqan3::format_bam{}, seqan3::fields<seqan3::field::alignment,
                                                                                      seqan3::field::tags>{}};
        ASSERT_THROW(fin.begin(), seqan3::format_error);
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

        seqan3::sam_file_input fin{stream, seqan3::format_bam{}};
        ASSERT_THROW(fin.begin(), seqan3::format_error);
    }
}

TEST_F(bam_format, too_long_cigar_string_write)
{
    // create an alignment resulting more than 65535 cigar elements
    // -------------------------------------------------------------------------
    auto read = seqan3::views::repeat_n('T'_dna5,  70'000);
    auto ref  = seqan3::views::repeat_n('A'_dna5, 2 * read.size() - 1);

    auto gapped_ref  = seqan3::gap_decorator{ref};

    // a gap_decorator on a repeat_n view also works but is slow when inserting gaps.
    std::vector<seqan3::gapped<seqan3::dna5>> gapped_read;
    gapped_read.reserve(2 * read.size());
    // create gap of length one every second character => T-T-T-T-T-T...
    for (auto chr : read)
        gapped_read.push_back(chr), gapped_read.push_back(seqan3::gap{});
    gapped_read.pop_back(); // remove last seqan3::gap

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
        (seqan3::views::repeat_n('\x88', (read.size() + 1) / 2) | seqan3::views::to<std::string> /*seq */) +
        (seqan3::views::repeat_n('\xFF',  read.size())          | seqan3::views::to<std::string> /*qual*/) +
        std::string /*the beginning*/
        {
            '\x43', '\x47', '\x5A' /*tag info: CGZ*/
        };

    for (size_t i = 0; i < read.size() - 1; ++i)
        expected.append("\x31\x4D\x31\x44"/*1M1D*/);
    expected.append("\x31\x4D"/*1M*/);
    expected.push_back('\x00' /*\0*/);

    std::ostringstream os{};

    seqan3::sam_file_header header{std::vector<std::string>{this->ref_id}};
    header.ref_id_info.push_back({ref.size(), ""});
    header.ref_dict[this->ref_id] = 0;

    {
        seqan3::sam_file_output fout{os, seqan3::format_bam{}, seqan3::fields<seqan3::field::header_ptr,
                                                                              seqan3::field::id,
                                                                              seqan3::field::seq,
                                                                              seqan3::field::ref_id,
                                                                              seqan3::field::ref_offset,
                                                                              seqan3::field::alignment,
                                                                              seqan3::field::mapq>{}};

        fout.emplace_back(&header, std::string{"long_read"}, read, 0, 0, alignment, 255);
    }

    os.flush();

    EXPECT_TRUE(os.str() == expected); // do not use EXPECT_EQ because if this fails the output will be huge :D
}

// https://github.com/seqan/seqan3/issues/2417
TEST_F(bam_format, issue2417)
{
    std::string const input{
        // @HD	VN:1.6
        // @SQ	SN:ref	LN:1904
        // read1	117	ref	1	0	*	=	1	0	ACGTA	IIIII
        '\x42', '\x41', '\x4D', '\x01', '\x1E', '\x00', '\x00', '\x00', '\x40', '\x48', '\x44', '\x09', '\x56', '\x4E',
        '\x3A', '\x31', '\x2E', '\x36', '\x0A', '\x40', '\x53', '\x51', '\x09', '\x53', '\x4E', '\x3A', '\x72', '\x65',
        '\x66', '\x09', '\x4C', '\x4E', '\x3A', '\x31', '\x39', '\x30', '\x34', '\x0A', '\x01', '\x00', '\x00', '\x00',
        '\x04', '\x00', '\x00', '\x00', '\x72', '\x65', '\x66', '\x00', '\x70', '\x07', '\x00', '\x00', '\x2E', '\x00',
        '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x06', '\x00', '\x49', '\x12',
        '\x00', '\x00', '\x75', '\x00', '\x05', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00',
        '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x72', '\x65', '\x61', '\x64', '\x31', '\x00', '\x12', '\x48',
        '\x10', '\x28', '\x28', '\x28', '\x28', '\x28'
    };

    std::istringstream stream{input};

    seqan3::sam_file_input fin{stream, seqan3::format_bam{}, seqan3::fields<seqan3::field::id,
                                                                            seqan3::field::alignment>{}};

    std::vector<seqan3::gapped<seqan3::dna5>> const empty_sequence{};

    size_t num_records{0u};

    // In 2417, the sequence was not consumed. Thus, wrong bytes were read for the following records.
    // With the chosen `input` this also means that there will be more than 1 record in the alignment file.
    // Hence, we need the for loop even though there is only 1 record.
    for (auto && [id, alignment] : fin)
    {
        ++num_records;
        EXPECT_RANGE_EQ(id, std::string{"read1"});
        EXPECT_RANGE_EQ(std::get<0>(alignment), empty_sequence);
        EXPECT_RANGE_EQ(std::get<1>(alignment), empty_sequence);
    }

    EXPECT_EQ(num_records, 1u);
}
