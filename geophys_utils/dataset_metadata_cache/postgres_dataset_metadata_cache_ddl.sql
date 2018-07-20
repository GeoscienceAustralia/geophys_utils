--
-- PostgreSQL database dump
--

-- Dumped from database version 10.1
-- Dumped by pg_dump version 10.4

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET client_min_messages = warning;
SET row_security = off;

--
-- Name: plpgsql; Type: EXTENSION; Schema: -; Owner: 
--

CREATE EXTENSION IF NOT EXISTS plpgsql WITH SCHEMA pg_catalog;


--
-- Name: EXTENSION plpgsql; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION plpgsql IS 'PL/pgSQL procedural language';


SET default_tablespace = '';

SET default_with_oids = false;

--
-- Name: dataset; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.dataset (
    dataset_id bigint NOT NULL,
    dataset_title character varying(128) NOT NULL,
    survey_id bigint,
    longitude_min double precision NOT NULL,
    longitude_max double precision NOT NULL,
    latitude_min double precision NOT NULL,
    latitude_max double precision NOT NULL,
    convex_hull_polygon text,
    metadata_uuid character(36)
);


ALTER TABLE public.dataset OWNER TO postgres;

--
-- Name: TABLE dataset; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.dataset IS 'Table of datasets to search';


--
-- Name: COLUMN dataset.dataset_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.dataset.dataset_id IS 'Numeric primary key for dataset';


--
-- Name: COLUMN dataset.dataset_title; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.dataset.dataset_title IS 'Dataset title string';


--
-- Name: COLUMN dataset.survey_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.dataset.survey_id IS 'Foreign key to survey (Null allowed)';


--
-- Name: COLUMN dataset.longitude_min; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.dataset.longitude_min IS 'Minimum longitude';


--
-- Name: COLUMN dataset.longitude_max; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.dataset.longitude_max IS 'Maximum longitude';


--
-- Name: COLUMN dataset.latitude_min; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.dataset.latitude_min IS 'Minimum latitude';


--
-- Name: COLUMN dataset.latitude_max; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.dataset.latitude_max IS 'Maximum latitude';


--
-- Name: COLUMN dataset.convex_hull_polygon; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.dataset.convex_hull_polygon IS 'Definition of convex hull polygon for dataset';


--
-- Name: dataset_dataset_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.dataset_dataset_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.dataset_dataset_id_seq OWNER TO postgres;

--
-- Name: dataset_dataset_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.dataset_dataset_id_seq OWNED BY public.dataset.dataset_id;


--
-- Name: dataset_keyword; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.dataset_keyword (
    dataset_id bigint NOT NULL,
    keyword_id bigint NOT NULL
);


ALTER TABLE public.dataset_keyword OWNER TO postgres;

--
-- Name: TABLE dataset_keyword; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.dataset_keyword IS 'Table expressing many-many relationships between datasets and keywords';


--
-- Name: COLUMN dataset_keyword.dataset_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.dataset_keyword.dataset_id IS 'Foreign key to dataset';


--
-- Name: COLUMN dataset_keyword.keyword_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.dataset_keyword.keyword_id IS 'Foreign key to keyword';


--
-- Name: distribution; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.distribution (
    distribution_id bigint NOT NULL,
    dataset_id bigint NOT NULL,
    distribution_url character varying(256) NOT NULL,
    protocol_id bigint NOT NULL
);


ALTER TABLE public.distribution OWNER TO postgres;

--
-- Name: TABLE distribution; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.distribution IS 'Table of online distributions';


--
-- Name: COLUMN distribution.distribution_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.distribution.distribution_id IS 'Numeric primary key for distribution';


--
-- Name: COLUMN distribution.dataset_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.distribution.dataset_id IS 'Foreign key to dataset';


--
-- Name: COLUMN distribution.distribution_url; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.distribution.distribution_url IS 'URL of distribution';


--
-- Name: COLUMN distribution.protocol_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.distribution.protocol_id IS 'Foreign key to protocol';


--
-- Name: distribution_distribution_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.distribution_distribution_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.distribution_distribution_id_seq OWNER TO postgres;

--
-- Name: distribution_distribution_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.distribution_distribution_id_seq OWNED BY public.distribution.distribution_id;


--
-- Name: keyword_keyword_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.keyword_keyword_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.keyword_keyword_id_seq OWNER TO postgres;

--
-- Name: keyword; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.keyword (
    keyword_id bigint DEFAULT nextval('public.keyword_keyword_id_seq'::regclass) NOT NULL,
    keyword_value character varying(64) NOT NULL,
    keyword_url character varying(256)
);


ALTER TABLE public.keyword OWNER TO postgres;

--
-- Name: TABLE keyword; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.keyword IS 'Table of keywords for datasets';


--
-- Name: COLUMN keyword.keyword_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.keyword.keyword_id IS 'Numeric primary key for keyword';


--
-- Name: COLUMN keyword.keyword_value; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.keyword.keyword_value IS 'Keyword string';


--
-- Name: COLUMN keyword.keyword_url; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.keyword.keyword_url IS 'URL of keyword (if known)';


--
-- Name: protocol; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.protocol (
    protocol_id bigint NOT NULL,
    protocol_value character varying(32) NOT NULL
);


ALTER TABLE public.protocol OWNER TO postgres;

--
-- Name: COLUMN protocol.protocol_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.protocol.protocol_id IS 'Numeric primary key for distribution protocol';


--
-- Name: COLUMN protocol.protocol_value; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.protocol.protocol_value IS 'Protocol value string';


--
-- Name: protocol_protocol_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.protocol_protocol_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.protocol_protocol_id_seq OWNER TO postgres;

--
-- Name: protocol_protocol_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.protocol_protocol_id_seq OWNED BY public.protocol.protocol_id;


--
-- Name: survey; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.survey (
    survey_id bigint NOT NULL,
    ga_survey_id character varying(16) NOT NULL,
    survey_name character varying(128)
);


ALTER TABLE public.survey OWNER TO postgres;

--
-- Name: TABLE survey; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.survey IS 'Table of survey metadata';


--
-- Name: COLUMN survey.survey_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.survey.survey_id IS 'Numeric primary key for survey entity';


--
-- Name: COLUMN survey.ga_survey_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.survey.ga_survey_id IS 'GA survey ID string';


--
-- Name: COLUMN survey.survey_name; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.survey.survey_name IS 'Survey name string';


--
-- Name: survey_survey_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.survey_survey_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.survey_survey_id_seq OWNER TO postgres;

--
-- Name: survey_survey_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.survey_survey_id_seq OWNED BY public.survey.survey_id;


--
-- Name: dataset dataset_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.dataset ALTER COLUMN dataset_id SET DEFAULT nextval('public.dataset_dataset_id_seq'::regclass);


--
-- Name: distribution distribution_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.distribution ALTER COLUMN distribution_id SET DEFAULT nextval('public.distribution_distribution_id_seq'::regclass);


--
-- Name: protocol protocol_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.protocol ALTER COLUMN protocol_id SET DEFAULT nextval('public.protocol_protocol_id_seq'::regclass);


--
-- Name: survey survey_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.survey ALTER COLUMN survey_id SET DEFAULT nextval('public.survey_survey_id_seq'::regclass);


--
-- Name: dataset_keyword dataset_keyword_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.dataset_keyword
    ADD CONSTRAINT dataset_keyword_pkey PRIMARY KEY (dataset_id, keyword_id);


--
-- Name: dataset dataset_metadata_uuid_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.dataset
    ADD CONSTRAINT dataset_metadata_uuid_key UNIQUE (metadata_uuid);


--
-- Name: dataset dataset_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.dataset
    ADD CONSTRAINT dataset_pkey PRIMARY KEY (dataset_id);


--
-- Name: distribution distribution_distribution_url_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.distribution
    ADD CONSTRAINT distribution_distribution_url_key UNIQUE (distribution_url);


--
-- Name: distribution distribution_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.distribution
    ADD CONSTRAINT distribution_pkey PRIMARY KEY (dataset_id, distribution_id);


--
-- Name: keyword keyword_keyword_value_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.keyword
    ADD CONSTRAINT keyword_keyword_value_key UNIQUE (keyword_value);


--
-- Name: keyword keyword_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.keyword
    ADD CONSTRAINT keyword_pkey PRIMARY KEY (keyword_id);


--
-- Name: dataset latitude_max_constraint; Type: CHECK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE public.dataset
    ADD CONSTRAINT latitude_max_constraint CHECK (((('-90'::integer)::double precision <= latitude_max) AND (latitude_max <= (90)::double precision))) NOT VALID;


--
-- Name: CONSTRAINT latitude_max_constraint ON dataset; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON CONSTRAINT latitude_max_constraint ON public.dataset IS 'Range check constraint for maximum latitude';


--
-- Name: dataset latitude_min_constraint; Type: CHECK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE public.dataset
    ADD CONSTRAINT latitude_min_constraint CHECK (((('-90'::integer)::double precision <= latitude_min) AND (latitude_min <= (90)::double precision))) NOT VALID;


--
-- Name: CONSTRAINT latitude_min_constraint ON dataset; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON CONSTRAINT latitude_min_constraint ON public.dataset IS 'Range check constraint for minimum latitude';


--
-- Name: dataset longitude_max_constraint; Type: CHECK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE public.dataset
    ADD CONSTRAINT longitude_max_constraint CHECK (((('-180'::integer)::double precision < longitude_max) AND (longitude_max <= (180)::double precision))) NOT VALID;


--
-- Name: CONSTRAINT longitude_max_constraint ON dataset; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON CONSTRAINT longitude_max_constraint ON public.dataset IS 'Range check constraint for maximum longitude';


--
-- Name: dataset longitude_min_constraint; Type: CHECK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE public.dataset
    ADD CONSTRAINT longitude_min_constraint CHECK (((('-180'::integer)::double precision < longitude_min) AND (longitude_min <= (180)::double precision))) NOT VALID;


--
-- Name: CONSTRAINT longitude_min_constraint ON dataset; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON CONSTRAINT longitude_min_constraint ON public.dataset IS 'Range check constraint for minimum longitude';


--
-- Name: protocol protocol_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.protocol
    ADD CONSTRAINT protocol_pkey PRIMARY KEY (protocol_id);


--
-- Name: protocol protocol_protocol_value_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.protocol
    ADD CONSTRAINT protocol_protocol_value_key UNIQUE (protocol_value);


--
-- Name: survey survey_ga_survey_id_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.survey
    ADD CONSTRAINT survey_ga_survey_id_key UNIQUE (ga_survey_id);


--
-- Name: survey survey_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.survey
    ADD CONSTRAINT survey_pkey PRIMARY KEY (survey_id);


--
-- Name: dataset_latitude_max_idx; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX dataset_latitude_max_idx ON public.dataset USING btree (latitude_max);


--
-- Name: dataset_latitude_min_idx; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX dataset_latitude_min_idx ON public.dataset USING btree (latitude_min);


--
-- Name: dataset_longitude_max_idx; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX dataset_longitude_max_idx ON public.dataset USING btree (longitude_max);


--
-- Name: dataset_longitude_min_idx; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX dataset_longitude_min_idx ON public.dataset USING btree (longitude_min);


--
-- Name: fki_dataset_keyword_dataset_id; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX fki_dataset_keyword_dataset_id ON public.dataset_keyword USING btree (dataset_id);


--
-- Name: fki_dataset_keyword_keyword_id; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX fki_dataset_keyword_keyword_id ON public.dataset_keyword USING btree (keyword_id);


--
-- Name: fki_dataset_survey_id; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX fki_dataset_survey_id ON public.dataset USING btree (survey_id);


--
-- Name: fki_distribution_dataset_id; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX fki_distribution_dataset_id ON public.distribution USING btree (dataset_id);


--
-- Name: fki_distribution_protocol_id; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX fki_distribution_protocol_id ON public.distribution USING btree (protocol_id);


--
-- Name: keyword_keyword_value_idx; Type: INDEX; Schema: public; Owner: postgres
--

CREATE UNIQUE INDEX keyword_keyword_value_idx ON public.keyword USING btree (keyword_value);


--
-- Name: dataset_keyword dataset_keyword_dataset_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.dataset_keyword
    ADD CONSTRAINT dataset_keyword_dataset_id_fkey FOREIGN KEY (dataset_id) REFERENCES public.dataset(dataset_id) ON UPDATE CASCADE;


--
-- Name: dataset_keyword dataset_keyword_keyword_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.dataset_keyword
    ADD CONSTRAINT dataset_keyword_keyword_id_fkey FOREIGN KEY (keyword_id) REFERENCES public.keyword(keyword_id) ON UPDATE CASCADE;


--
-- Name: dataset dataset_survey_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.dataset
    ADD CONSTRAINT dataset_survey_id_fkey FOREIGN KEY (survey_id) REFERENCES public.survey(survey_id) ON UPDATE CASCADE;


--
-- Name: distribution distribution_dataset_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.distribution
    ADD CONSTRAINT distribution_dataset_id_fkey FOREIGN KEY (dataset_id) REFERENCES public.dataset(dataset_id) ON UPDATE CASCADE;


--
-- Name: distribution distribution_protocol_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.distribution
    ADD CONSTRAINT distribution_protocol_id_fkey FOREIGN KEY (protocol_id) REFERENCES public.protocol(protocol_id) ON UPDATE CASCADE;


--
-- Name: TABLE dataset; Type: ACL; Schema: public; Owner: postgres
--

GRANT ALL ON TABLE public.dataset TO db_users;


--
-- Name: SEQUENCE dataset_dataset_id_seq; Type: ACL; Schema: public; Owner: postgres
--

GRANT ALL ON SEQUENCE public.dataset_dataset_id_seq TO db_users;


--
-- Name: TABLE dataset_keyword; Type: ACL; Schema: public; Owner: postgres
--

GRANT ALL ON TABLE public.dataset_keyword TO db_users;


--
-- Name: TABLE distribution; Type: ACL; Schema: public; Owner: postgres
--

GRANT ALL ON TABLE public.distribution TO db_users;


--
-- Name: SEQUENCE distribution_distribution_id_seq; Type: ACL; Schema: public; Owner: postgres
--

GRANT ALL ON SEQUENCE public.distribution_distribution_id_seq TO db_users;


--
-- Name: SEQUENCE keyword_keyword_id_seq; Type: ACL; Schema: public; Owner: postgres
--

GRANT ALL ON SEQUENCE public.keyword_keyword_id_seq TO db_users;


--
-- Name: TABLE keyword; Type: ACL; Schema: public; Owner: postgres
--

GRANT ALL ON TABLE public.keyword TO db_users;


--
-- Name: TABLE protocol; Type: ACL; Schema: public; Owner: postgres
--

GRANT ALL ON TABLE public.protocol TO db_users;


--
-- Name: SEQUENCE protocol_protocol_id_seq; Type: ACL; Schema: public; Owner: postgres
--

GRANT ALL ON SEQUENCE public.protocol_protocol_id_seq TO db_users;


--
-- Name: TABLE survey; Type: ACL; Schema: public; Owner: postgres
--

GRANT ALL ON TABLE public.survey TO db_users;


--
-- Name: SEQUENCE survey_survey_id_seq; Type: ACL; Schema: public; Owner: postgres
--

GRANT ALL ON SEQUENCE public.survey_survey_id_seq TO db_users;


--
-- PostgreSQL database dump complete
--

