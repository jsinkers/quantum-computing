module ComplexNumbers
  (
    ComplexNumber,
    Re,
    Im,
    Rad,
    Theta,
    ComplexVector,
    ComplexMatrix,
    complexAdd,
    complexMult,
    complexSub,
    complexModulus,
    complexConjugate,
    complexInverse,
    complexPower,
    complexRoots,
    complexDiv,
    complexRe,
    complexIm,
    complexRadius,
    complexPhase,
    toPolar,
    toCartesian,
    cvAdd,
    cvScalarMult,
    cvInverse,
    cvInnerProd,
    cvDot,
    cvConjugate,
    cvNorm,
    cvDist,
    cmIdentity,
    cmAdd,
    cmScalarMult,
    cmInverse,
    cmTranspose,
    cmDagger,
    cmAdjoint,
    cmConjugate,
    cmMult,
    cmInnerProd,
    cmTrace,
    cmNorm,
    printComplexMatrix,
    showComplexMatrix
   ) where 

import Data.List (transpose, intercalate)

-- todo: introduce a type class for inner space product etc for polymorphic
-- functions, instead of functions particular to vectors and matrices

-- data types
-- ============================================================================
type Re = Double
type Im = Double

type Rad = Double
type Theta = Double

data ComplexNumber = Cartesian Re Im | Euler Rad Theta


type ComplexVector = [ComplexNumber]

-- m x n complex matrix: contains m lists of n elements
type ComplexMatrix = [[ComplexNumber]]

instance Show ComplexNumber where
    show (Cartesian r i) = show r ++ " + " ++ show i ++ "i"
    show (Euler r t) = show r ++ "e^{i" ++ show t ++ "}"
instance Eq ComplexNumber where
    a == b = complexRe a' == complexRe b' && complexIm a' == complexIm b'
        where a' = toCartesian a
              b' = toCartesian b

-- complex number algebra
-- ============================================================================
-- complex zero number
complexZero :: ComplexNumber
complexZero = Cartesian 0 0

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

-- determine inverse of a complex number (i.e. z^-1 = 1/z) - not to be confused with additive inverse
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

-- complex vector space operations
-- ============================================================================
cvAdd :: ComplexVector -> ComplexVector -> ComplexVector
cvAdd cv1s cv2s = if length cv1s == length cv2s 
                     then [complexAdd c1 c2 | (c1, c2) <- zip cv1s cv2s]
                     else error "ComplexVector size mismatch"

cvScalarMult :: ComplexNumber -> ComplexVector -> ComplexVector
cvScalarMult s = map (complexMult s) 

cvInverse :: ComplexVector -> ComplexVector
cvInverse = map (complexSub (Euler 0 0)) 

-- helper for matrix calculations
cvProduct :: ComplexVector -> ComplexVector -> ComplexNumber
cvProduct cv1 cv2 =  foldl complexAdd complexZero $ zipWith complexMult cv1 cv2

-- inner product/dot product
cvInnerProd :: ComplexVector -> ComplexVector -> ComplexNumber
cvInnerProd v1 v2 = if length v1 == length v2
                       then cvProduct (cvConjugate v1) v2
                       else error "ComplexVector size mismatch"

-- dot product: alias for inner product
cvDot :: ComplexVector -> ComplexVector -> ComplexNumber
cvDot = cvInnerProd

-- take conjugate of a complex vector
cvConjugate :: ComplexVector -> ComplexVector
cvConjugate = map complexConjugate 

-- determine norm of a complex vector
cvNorm :: ComplexVector -> Double
cvNorm v = sqrt $ complexRe $ cvInnerProd v v 

-- determine distance between two complex vectors
cvDist :: ComplexVector -> ComplexVector -> Double
cvDist v1 v2 = cvNorm $ cvAdd v1 $ cvInverse v2

-- determine if two complex vectors are orthogonal
cvOrthogonal :: ComplexVector -> ComplexVector -> Bool
cvOrthogonal v1 v2 = cvInnerProd v1 v2 == complexZero

-- complex matrix functions
-- ============================================================================
-- complex identity matrix I_n
cmIdentity :: Int -> ComplexMatrix
cmIdentity n = [[if i == j then Cartesian 1 0 else complexZero | j <- [0..n-1]] | i <- [0..n-1]]

-- determine result of addition of 2 complex matrices
cmAdd :: ComplexMatrix -> ComplexMatrix -> ComplexMatrix
cmAdd cm1s cm2s = if length cm1s == length cm2s 
                     then [cvAdd c1 c2 | (c1, c2) <- zip cm1s cm2s]
                     else error "ComplexMatrix size mismatch"

-- determine scalar multiple of a complex matrix
cmScalarMult :: ComplexNumber -> ComplexMatrix -> ComplexMatrix
cmScalarMult s = map (cvScalarMult s)

-- determine additive inverse of a complex matrix
cmInverse :: ComplexMatrix -> ComplexMatrix
cmInverse = map cvInverse

-- determine transpose of a complex matrix
cmTranspose :: ComplexMatrix -> ComplexMatrix
cmTranspose = transpose

-- determine conjugate of a complex matrix
cmConjugate :: ComplexMatrix -> ComplexMatrix
cmConjugate = map (map complexConjugate) 

-- determine adjoint/dagger
cmDagger :: ComplexMatrix -> ComplexMatrix
cmDagger = cmConjugate . cmTranspose

cmAdjoint :: ComplexMatrix -> ComplexMatrix
cmAdjoint = cmDagger

showComplexMatrix :: ComplexMatrix -> String
showComplexMatrix cm = intercalate "\n" $ map show cm 

printComplexMatrix :: ComplexMatrix -> IO ()
printComplexMatrix cm = do putStr $ showComplexMatrix cm
                           putStrLn ""

-- multiply 2 complex matrices
-- note no size checking is performed
cmMult :: ComplexMatrix -> ComplexMatrix -> ComplexMatrix
cmMult cm1 cm2 = [[cvProduct cv1 cv2 | cv2 <- cm2'] | cv1 <- cm1]
    where cm2' = cmTranspose cm2

-- Action of application of complex matrix A to complex vector v
cmAction :: ComplexMatrix -> ComplexVector -> ComplexVector
cmAction a v = [cvProduct ra v | ra <- a]

-- determine inner product of 2 complex matrices
cmInnerProd :: ComplexMatrix -> ComplexMatrix -> ComplexNumber
cmInnerProd m1 m2 = cmTrace $ cmMult (cmDagger m1) m2

-- determine trace of a complex matrix
cmTrace :: ComplexMatrix -> ComplexNumber
cmTrace cm = foldl complexAdd complexZero $ [r !! n | (r, n) <- zip cm [0..]]

-- determine norm of a complex matrix
cmNorm :: ComplexMatrix -> Double
cmNorm m = sqrt $ complexRe $ cmInnerProd m m

-- misc functions
kroneckerDelta :: Int -> Int -> Int
kroneckerDelta j k
    | j == k    = 1
    | otherwise = 0

-- determine if a matrix is hermitian
cmHermitian :: ComplexMatrix -> Bool
cmHermitian m = if cmSquare m 
                   then cmDagger m == m
                   else error "ComplexMatrix is not square"

-- determine if a complex matrix is square
cmSquare :: ComplexMatrix -> Bool
cmSquare m = all (== length m) $ map length m

-- determine if a complex matrix is unitary
cmUnitary :: ComplexMatrix -> Bool
cmUnitary m = if cmSquare m
                 then (cmMult m $ cmDagger m) == (cmIdentity $ length m)
                 else error "ComplexMatrix is not square"

-- determine the tensor product of two matrices
cmTensorProd :: ComplexMatrix -> ComplexMatrix -> ComplexMatrix
cmTensorProd cm1 cm2 = [[complexMult x y | x <- p, y <- q] | p <- cm1, q <- cm2]

-- test matrix
cm1 :: ComplexMatrix
cm1 = [[Cartesian 6 (-3), Cartesian 2 12,  Cartesian 0 (-19) ],
       [Cartesian 0 0,    Cartesian 5 2.1, Cartesian 17 0    ],
       [Cartesian 1 0,    Cartesian 2 5,   Cartesian 3 (-4.5)]]

cm2 :: ComplexMatrix
cm2 = [[Cartesian 1 0, Cartesian 0 0,  Cartesian 0 0 ],
                 [Cartesian 0 0,    Cartesian 1 0.0, Cartesian 0 0    ],
                 [Cartesian 0 0,    Cartesian 0 0,   Cartesian 1 0]]

cm3 :: ComplexMatrix
cm3 = [[Cartesian 1 0, Cartesian 1 1],
       [Cartesian 1 1, Cartesian 1 0]]

cm4 :: ComplexMatrix
cm4 = [[Cartesian 2 1, Cartesian 4 1],
       [Cartesian 3 0, Cartesian 5 0]]

cm5 :: ComplexMatrix
cm5 = [[Cartesian 3 2, Cartesian 5 (-1),  Cartesian 0 2 ],
       [Cartesian 0 0,    Cartesian 12 0, Cartesian 6 (-3)    ],
       [Cartesian 2 0,    Cartesian 4 4,   Cartesian 9 3]]

cm6 :: ComplexMatrix
cm6 = [[Cartesian 1 0, Cartesian 3 4,  Cartesian 5 (-7) ],
       [Cartesian 10 2,    Cartesian 6 0, Cartesian 2 5    ],
       [Cartesian 0 0,    Cartesian 1 0,   Cartesian 2 9]]
