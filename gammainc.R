gammainc <- function (x,a) {
		if (x == 0 || a == 0 || abs(a)>170) {gin <- NA}
		if (x > 0 && abs(a) <= 170) {
			xam <- -x+a*log(x)
			if (is.na(xam) || abs(xam) > 700) {gin <- NA}
			gin <- 0
			if (x <= 1+a) {
				s <- 1/a
				r <- s
				for (k in 1:60) {
					r <- r * x/(a + k)
					s <- s + r
					if (abs(r/s) < 1e-15)
					break
				}
				gin <- exp(xam) * s
			} else {
				t0 <- 0
				for (k in 60:1) {
					t0 <- (k - a)/(1 + k/(x + t0))
				}
				gin <- gamma(a) - exp(xam)/(x + t0)
			}
			gin
		}
	}
