# Function to perform terrain correction on Sentinel-1 images
def terrain_correction(image):
    img_geom = image.geometry()
    srtm = ee.Image('USGS/SRTMGL1_003').clip(img_geom)  # 30m SRTM
    sigma0_pow = ee.Image.constant(10).pow(image.divide(10.0))

    # Radar geometry
    theta_i = image.select('angle')
    phi_i = ee.Terrain.aspect(theta_i) \
        .reduceRegion(ee.Reducer.mean(), theta_i.get('system:footprint'), 1000) \
        .get('aspect')

    # Terrain geometry
    alpha_s = ee.Terrain.slope(srtm).select('slope')
    phi_s = ee.Terrain.aspect(srtm).select('aspect')

    # Model geometry
    phi_r = ee.Image.constant(phi_i).subtract(phi_s)

    # Convert all to radians
    phi_r_rad = phi_r.multiply((math.pi) / 180)
    alpha_s_rad = alpha_s.multiply((math.pi) / 180)
    theta_i_rad = theta_i.multiply((math.pi) / 180)
    ninety_rad = ee.Image.constant(90).multiply((math.pi) / 180)

    # Slope steepness in range
    alpha_r = alpha_s_rad.tan().multiply(phi_r_rad.cos()).atan()

    # Slope steepness in azimuth
    alpha_az = alpha_s_rad.tan().multiply(phi_r_rad.sin()).atan()

    # Local incidence angle
    theta_lia = alpha_az.cos().multiply(theta_i_rad.subtract(alpha_r).cos()).acos()
    theta_lia_deg = theta_lia.multiply(180 / math.pi)

    # Gamma nought flat
    gamma0 = sigma0_pow.divide(theta_i_rad.cos())
    gamma0_db = ee.Image.constant(10).multiply(gamma0.log10())
    ratio_1 = gamma0_db.select('VV').subtract(gamma0_db.select('VH'))

    # Volumetric Model
    nominator = ninety_rad.subtract(theta_i_rad).add(alpha_r).tan()
    denominator = ninety_rad.subtract(theta_i_rad).tan()
    vol_model = nominator.divide(denominator).abs()

    # Apply model
    gamma0_volume = gamma0.divide(vol_model)
    gamma0_volume_db = ee.Image.constant(10).multiply(gamma0_volume.log10())

    # Layover/shadow mask
    alpha_r_deg = alpha_r.multiply(180 / math.pi)
    layover = alpha_r_deg.lt(theta_i)

    # Shadow where LIA > 90
    shadow = theta_lia_deg.lt(85)

    # Ratio for RGB visualization
    ratio = gamma0_volume_db.select('VV').subtract(gamma0_volume_db.select('VH'))

    output = gamma0_volume_db.addBands(ratio).addBands(alpha_r).addBands(phi_s).addBands(theta_i_rad) \
        .addBands(layover).addBands(shadow).addBands(gamma0_db).addBands(ratio_1)

    return image.addBands(
        output.select(['VV', 'VH'], ['VV', 'VH']),
        None,
        True
    )
