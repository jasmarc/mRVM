require "rubygems"
require "Matrix"
require "rbgsl"

class Matrix
  def []=(i,j,x)
    @rows[i][j]=x
  end
end

y = Matrix[[1, 8, 2, 6],
           [4, 5, 1, 0],
           [3, 3, 4, 2]]

w = Matrix[[3, 1, 4, 0],
           [4, 3, 9, 2],
           [3, 1, 5, 7]]

k = Matrix[[1, 2, 9],
           [4, 6, 5],
           [3, 7, 4]]

classes = Vector[2, 1, 0]

N = y.row_size
C = y.column_size

0.upto(N-1) do |n|
  puts "=" * 20
  i = classes[n]
  puts "i = #{i}"
  puts "n = #{n}"
  k_n = k.row(n)
  puts "k_n = #{k_n}"
  0.upto(C-1) do |c|
    w_c = w.column(c)
    w_i = w.column(classes[n])
    wckn = w_c.inner_product k_n
    wikn = w_i.inner_product k_n
    puts "\t" + "="*20
    puts "\tc = #{c}"
    puts "\tw_c = #{w_c}"
    puts "\tw_i = #{w_i}"
    puts "\tw_c*k_n = #{wckn}"
    puts "\tw_i*k_n = #{wikn}"
    puts "\ty_cn current value = #{y[n, c]}"
    if c != i
      r = GSL::Rng.alloc(GSL::Rng::TAUS, 1)
      numerator   = 0
      denominator = 0
      0.upto(1000) do
        u = r.gaussian
        num = GSL::Ran::gaussian_pdf(wckn - wikn)
        den = GSL::Cdf::gaussian_P(u + wikn - wckn)
        0.upto(C-1) do |j|
          if(j != i && j != c)
            w_j = w.column(j)
            wjkn = w_j.inner_product k_n
            num *= GSL::Cdf::gaussian_P(u + wikn - wjkn)
            den *= GSL::Cdf::gaussian_P(u + wikn - wjkn)
          end
        end
        numerator   += num
        denominator += den
      end
      puts "\tnumerator = #{numerator}"
      puts "\tdenominator = #{denominator}"
      puts "\tvalue = #{numerator/denominator}"
      puts "\ty_cn = #{wckn} - #{numerator} / #{denominator}"
      puts "\ty_cn = #{wckn - (numerator / denominator)}"
      y[n, c] = wckn - (numerator / denominator)
    else
      val = wckn
      0.upto(C-1) do |j|
        if(j != i)
          w_j = w.column(j)
          wjkn = w_j.inner_product k_n
          puts "\twjkn = #{wjkn}"
          puts "\ty[n, j] = #{y[n, j]}"
          puts "\tval = #{val}"
          val -= (y[n, j] - wjkn)
        end
      end
      puts "\tval = #{val}"
      y[n, c] = val
    end
  end
end

# r = GSL::Rng.alloc(GSL::Rng::TAUS, 1)
# 0.upto(10) do
#   puts r.gaussian
# end

# -4.upto(4) do |x|
#   puts GSL::Ran::gaussian_pdf(x)
# end

# -5.upto(5) do |x|
#   puts "CDF #{GSL::Cdf::gaussian_P(x)}"
# end