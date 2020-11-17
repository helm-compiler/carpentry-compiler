// SExp - A S-Expression Parser for C++
// Copyright (C) 2018 Ingo Ruhnke <grumbel@gmail.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <gtest/gtest.h>

#include <string>
#include <sstream>
#include <vector>

#include "float.hpp"

TEST(FloatTest, float2string)
{
  std::vector<std::tuple<float, std::string> > datas =
    {
      { 1.0f, "1" },
      { 1234.567f, "1234.567" },
      { 1.0e3f, "1000" },
      { 1234567e-4f, "123.4567" },
      { 1.234e13f, "1.234e+13" },
      { 1234.0e13f, "1.234e+16" },
    };

  for(const auto& data: datas)
  {
    std::ostringstream os;
    sexp::float2string(os, std::get<0>(data));
    ASSERT_EQ(os.str(), std::get<1>(data));
  }
}

TEST(FloatTest, string2float)
{
  ASSERT_EQ(1234.5678f, sexp::string2float("1234.5678"));
  ASSERT_EQ(1234.5678f, sexp::string2float("+1234.5678"));
  ASSERT_EQ(-1234.5678f, sexp::string2float("-1234.5678"));

  ASSERT_EQ(5.0f, sexp::string2float("5.0ef"));
  ASSERT_EQ(1234.5f, sexp::string2float("1.2345e3f"));
  ASSERT_EQ(12.345f, sexp::string2float("12345e-3f"));
}

/* EOF */
