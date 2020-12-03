-- complex number algebra
type Re = Double
type Im = Double

type Rad = Double
type Theta = Double

data ComplexNumber = Cartesian Re Im | Euler Rad Theta

instance Show ComplexNumber where
    show (Cartesian r i) = show r ++ " + " ++ show i ++ "i"
    show (Euler r t) = show r ++ "e^{i" ++ show t ++ "}"


complexAdd :: ComplexNumber -> ComplexNumber -> ComplexNumber
complexAdd c1 c2 = Cartesian (r1+r2) (i1+i2)
    where (Cartesian r1 i1) = toCartesian c1 
          (Cartesian r2 i2) = toCartesian c2

complexMult :: ComplexNumber -> ComplexNumber -> ComplexNumber
complexMult (Euler r1 t1) (Euler r2 t2) = Euler (r1*r2) (t1+t2)
complexMult c1 c2 = Cartesian (r1*r2 - i1*i2) (r1*i2 + r2*i1)
    where (Cartesian r1 i1) = toCartesian c1
          (Cartesian r2 i2) = toCartesian c2

complexSub :: ComplexNumber -> ComplexNumber -> ComplexNumber
complexSub c1 c2 = Cartesian (r1-r2) (i1-i2)
    where (Cartesian r1 i1) = toCartesian c1
          (Cartesian r2 i2) = toCartesian c2

complexModulus :: ComplexNumber -> Double
complexModulus (Cartesian r i) = sqrt (r^2 + i^2)
complexModulus (Euler r _) = r

-- determine complex conjugate of a complex number
complexConjugate :: ComplexNumber -> ComplexNumber
complexConjugate (Cartesian r1 i1) = Cartesian r1 (-i1)
complexConjugate (Euler r t) = Euler r (-t)

-- determine inverse of a complex number
complexInverse :: ComplexNumber -> ComplexNumber
complexInverse c = complexMult (Cartesian (1/modCsq) 0) (complexConjugate c)
    where modCsq = complexModulus c ^ 2

-- division of complex numbers
complexDiv :: ComplexNumber -> ComplexNumber -> ComplexNumber
complexDiv c1 c2 = complexMult c1 $ complexInverse c2

-- convert to polar form
toPolar :: ComplexNumber -> ComplexNumber
toPolar p@(Euler _ _) = p
toPolar c@(Cartesian r i) = Euler (complexModulus c) (atan2 i r)

-- convert to Cartesian form
toCartesian :: ComplexNumber -> ComplexNumber
toCartesian c@(Cartesian _ _) = c
toCartesian (Euler r theta) = Cartesian (r*cos theta) (r*sin theta)

-- getter for real part of a complex number
complexRe :: ComplexNumber -> Re
complexRe c = r
    where (Cartesian r _) = toCartesian c

-- getter for imaginary part of a complex number
complexIm :: ComplexNumber -> Im
complexIm c = i
    where (Cartesian _ i) = toCartesian c

-- getter for magnitude of a complex number
complexRadius :: ComplexNumber -> Rad
complexRadius = complexModulus

-- getter for phase of a complex number
complexPhase :: ComplexNumber -> Theta
complexPhase c = t
    where (Euler _ t) = toPolar c

complexPower :: Double -> ComplexNumber -> ComplexNumber
complexPower n c = Euler (r**n) (n*t)
    where (Euler r t) = toPolar c

-- need a function to convert angles to (-Pi, Pi] or [0,2Pi]


-- determine roots of unity
complexRoots :: Int -> [ComplexNumber]
complexRoots n = [Euler 1 (fromIntegral (2*k) * pi/fromIntegral n) | k <- [0..n-1]]
