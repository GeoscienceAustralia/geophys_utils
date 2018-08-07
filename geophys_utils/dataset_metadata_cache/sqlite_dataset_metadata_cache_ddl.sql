CREATE TABLE dataset (
    dataset_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL UNIQUE,
    dataset_title character varying(128) NOT NULL,
    survey_id INTEGER,
    longitude_min double precision NOT NULL CHECK((-180 < longitude_min) AND (longitude_min <= 180)),
    longitude_max double precision NOT NULL CHECK((-180 < longitude_max) AND (longitude_max <= 180)),
    latitude_min double precision NOT NULL CHECK((-90 <= latitude_min) AND (latitude_min <= 90)),
    latitude_max double precision NOT NULL CHECK((-90 <= latitude_max) AND (latitude_max <= 90)),
    convex_hull_polygon text,
    metadata_uuid character(36) NOT NULL UNIQUE,
    point_count INTEGER,
    FOREIGN KEY (survey_id) REFERENCES survey(survey_id) ON UPDATE CASCADE
);

CREATE TABLE dataset_keyword (
    dataset_id INTEGER NOT NULL,
    keyword_id INTEGER NOT NULL,
	FOREIGN KEY (dataset_id) REFERENCES dataset(dataset_id) ON UPDATE CASCADE,
	FOREIGN KEY (keyword_id) REFERENCES keyword(keyword_id) ON UPDATE CASCADE,
	UNIQUE (dataset_id, keyword_id)
);

CREATE TABLE distribution (
    distribution_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL UNIQUE,
    dataset_id INTEGER NOT NULL,
    distribution_url character varying(256) NOT NULL, -- Should be globally unique
    protocol_id INTEGER NOT NULL,
    FOREIGN KEY (dataset_id) REFERENCES dataset(dataset_id) ON UPDATE CASCADE,
    FOREIGN KEY (protocol_id) REFERENCES protocol(protocol_id) ON UPDATE CASCADE
	UNIQUE (dataset_id, distribution_url)
);

CREATE TABLE keyword (
    keyword_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL UNIQUE,
    keyword_value character varying(64) NOT NULL UNIQUE,
    keyword_url character varying(256)
);

CREATE TABLE protocol (
    protocol_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL UNIQUE,
    protocol_value character varying(32) NOT NULL UNIQUE
);

CREATE TABLE survey (
    survey_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL UNIQUE,
    ga_survey_id character varying(16) NOT NULL UNIQUE,
    survey_name character varying(128),
    start_date TIMESTAMP,
    end_date TIMESTAMP
);

CREATE INDEX fki_dataset_keyword_dataset_id ON dataset_keyword (dataset_id);

CREATE INDEX fki_dataset_keyword_keyword_id ON dataset_keyword (keyword_id);

CREATE INDEX fki_dataset_survey_id ON dataset (survey_id);

CREATE INDEX fki_distribution_dataset_id ON distribution (dataset_id);

CREATE INDEX fki_distribution_protocol_id ON distribution (protocol_id);



CREATE INDEX dataset_longitude_min_idx ON dataset (longitude_min);
CREATE INDEX dataset_longitude_max_idx ON dataset (longitude_max);
CREATE INDEX dataset_latitude_min_idx ON dataset (latitude_min);
CREATE INDEX dataset_latitude_max_idx ON dataset (latitude_max);

