--
-- PostgreSQL database dump
--

-- Dumped from database version 10.1
-- Dumped by pg_dump version 10.4

-- Started on 2018-07-20 08:55:38

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
-- TOC entry 1 (class 3079 OID 12924)
-- Name: plpgsql; Type: EXTENSION; Schema: -; Owner: 
--

CREATE EXTENSION IF NOT EXISTS plpgsql WITH SCHEMA pg_catalog;


--
-- TOC entry 2868 (class 0 OID 0)
-- Dependencies: 1
-- Name: EXTENSION plpgsql; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION plpgsql IS 'PL/pgSQL procedural language';


SET default_tablespace = '';

SET default_with_oids = false;

--
-- TOC entry 196 (class 1259 OID 29769)
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
-- TOC entry 2869 (class 0 OID 0)
-- Dependencies: 196
-- Name: TABLE dataset; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.dataset IS 'Table of datasets to search';


--
-- TOC entry 2870 (class 0 OID 0)
-- Dependencies: 196
-- Name: COLUMN dataset.dataset_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.dataset.dataset_id IS 'Numeric primary key for dataset';


--
-- TOC entry 2871 (class 0 OID 0)
-- Dependencies: 196
-- Name: COLUMN dataset.dataset_title; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.dataset.dataset_title IS 'Dataset title string';


--
-- TOC entry 2872 (class 0 OID 0)
-- Dependencies: 196
-- Name: COLUMN dataset.survey_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.dataset.survey_id IS 'Foreign key to survey (Null allowed)';


--
-- TOC entry 2873 (class 0 OID 0)
-- Dependencies: 196
-- Name: COLUMN dataset.longitude_min; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.dataset.longitude_min IS 'Minimum longitude';


--
-- TOC entry 2874 (class 0 OID 0)
-- Dependencies: 196
-- Name: COLUMN dataset.longitude_max; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.dataset.longitude_max IS 'Maximum longitude';


--
-- TOC entry 2875 (class 0 OID 0)
-- Dependencies: 196
-- Name: COLUMN dataset.latitude_min; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.dataset.latitude_min IS 'Minimum latitude';


--
-- TOC entry 2876 (class 0 OID 0)
-- Dependencies: 196
-- Name: COLUMN dataset.latitude_max; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.dataset.latitude_max IS 'Maximum latitude';


--
-- TOC entry 2877 (class 0 OID 0)
-- Dependencies: 196
-- Name: COLUMN dataset.convex_hull_polygon; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.dataset.convex_hull_polygon IS 'Definition of convex hull polygon for dataset';


--
-- TOC entry 197 (class 1259 OID 29775)
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
-- TOC entry 2878 (class 0 OID 0)
-- Dependencies: 197
-- Name: dataset_dataset_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.dataset_dataset_id_seq OWNED BY public.dataset.dataset_id;


--
-- TOC entry 198 (class 1259 OID 29777)
-- Name: dataset_keyword; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.dataset_keyword (
    dataset_id bigint NOT NULL,
    keyword_id bigint NOT NULL
);


ALTER TABLE public.dataset_keyword OWNER TO postgres;

--
-- TOC entry 2879 (class 0 OID 0)
-- Dependencies: 198
-- Name: TABLE dataset_keyword; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.dataset_keyword IS 'Table expressing many-many relationships between datasets and keywords';


--
-- TOC entry 2880 (class 0 OID 0)
-- Dependencies: 198
-- Name: COLUMN dataset_keyword.dataset_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.dataset_keyword.dataset_id IS 'Foreign key to dataset';


--
-- TOC entry 2881 (class 0 OID 0)
-- Dependencies: 198
-- Name: COLUMN dataset_keyword.keyword_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.dataset_keyword.keyword_id IS 'Foreign key to keyword';


--
-- TOC entry 199 (class 1259 OID 29780)
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
-- TOC entry 2882 (class 0 OID 0)
-- Dependencies: 199
-- Name: TABLE distribution; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.distribution IS 'Table of online distributions';


--
-- TOC entry 2883 (class 0 OID 0)
-- Dependencies: 199
-- Name: COLUMN distribution.distribution_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.distribution.distribution_id IS 'Numeric primary key for distribution';


--
-- TOC entry 2884 (class 0 OID 0)
-- Dependencies: 199
-- Name: COLUMN distribution.dataset_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.distribution.dataset_id IS 'Foreign key to dataset';


--
-- TOC entry 2885 (class 0 OID 0)
-- Dependencies: 199
-- Name: COLUMN distribution.distribution_url; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.distribution.distribution_url IS 'URL of distribution';


--
-- TOC entry 2886 (class 0 OID 0)
-- Dependencies: 199
-- Name: COLUMN distribution.protocol_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.distribution.protocol_id IS 'Foreign key to protocol';


--
-- TOC entry 200 (class 1259 OID 29783)
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
-- TOC entry 2887 (class 0 OID 0)
-- Dependencies: 200
-- Name: distribution_distribution_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.distribution_distribution_id_seq OWNED BY public.distribution.distribution_id;


--
-- TOC entry 201 (class 1259 OID 29785)
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
-- TOC entry 202 (class 1259 OID 29787)
-- Name: keyword; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.keyword (
    keyword_id bigint DEFAULT nextval('public.keyword_keyword_id_seq'::regclass) NOT NULL,
    keyword_value character varying(64) NOT NULL,
    keyword_url character varying(256)
);


ALTER TABLE public.keyword OWNER TO postgres;

--
-- TOC entry 2888 (class 0 OID 0)
-- Dependencies: 202
-- Name: TABLE keyword; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.keyword IS 'Table of keywords for datasets';


--
-- TOC entry 2889 (class 0 OID 0)
-- Dependencies: 202
-- Name: COLUMN keyword.keyword_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.keyword.keyword_id IS 'Numeric primary key for keyword';


--
-- TOC entry 2890 (class 0 OID 0)
-- Dependencies: 202
-- Name: COLUMN keyword.keyword_value; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.keyword.keyword_value IS 'Keyword string';


--
-- TOC entry 2891 (class 0 OID 0)
-- Dependencies: 202
-- Name: COLUMN keyword.keyword_url; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.keyword.keyword_url IS 'URL of keyword (if known)';


--
-- TOC entry 203 (class 1259 OID 29791)
-- Name: protocol; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.protocol (
    protocol_id bigint NOT NULL,
    protocol_value character varying(32) NOT NULL
);


ALTER TABLE public.protocol OWNER TO postgres;

--
-- TOC entry 2892 (class 0 OID 0)
-- Dependencies: 203
-- Name: COLUMN protocol.protocol_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.protocol.protocol_id IS 'Numeric primary key for distribution protocol';


--
-- TOC entry 2893 (class 0 OID 0)
-- Dependencies: 203
-- Name: COLUMN protocol.protocol_value; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.protocol.protocol_value IS 'Protocol value string';


--
-- TOC entry 204 (class 1259 OID 29794)
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
-- TOC entry 2894 (class 0 OID 0)
-- Dependencies: 204
-- Name: protocol_protocol_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.protocol_protocol_id_seq OWNED BY public.protocol.protocol_id;


--
-- TOC entry 205 (class 1259 OID 29796)
-- Name: survey; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.survey (
    survey_id bigint NOT NULL,
    ga_survey_id character varying(16) NOT NULL,
    survey_name character varying(128)
);


ALTER TABLE public.survey OWNER TO postgres;

--
-- TOC entry 2895 (class 0 OID 0)
-- Dependencies: 205
-- Name: TABLE survey; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.survey IS 'Table of survey metadata';


--
-- TOC entry 2896 (class 0 OID 0)
-- Dependencies: 205
-- Name: COLUMN survey.survey_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.survey.survey_id IS 'Numeric primary key for survey entity';


--
-- TOC entry 2897 (class 0 OID 0)
-- Dependencies: 205
-- Name: COLUMN survey.ga_survey_id; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.survey.ga_survey_id IS 'GA survey ID string';


--
-- TOC entry 2898 (class 0 OID 0)
-- Dependencies: 205
-- Name: COLUMN survey.survey_name; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.survey.survey_name IS 'Survey name string';


--
-- TOC entry 206 (class 1259 OID 29799)
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
-- TOC entry 2899 (class 0 OID 0)
-- Dependencies: 206
-- Name: survey_survey_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.survey_survey_id_seq OWNED BY public.survey.survey_id;


--
-- TOC entry 2699 (class 2604 OID 29801)
-- Name: dataset dataset_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.dataset ALTER COLUMN dataset_id SET DEFAULT nextval('public.dataset_dataset_id_seq'::regclass);


--
-- TOC entry 2704 (class 2604 OID 29802)
-- Name: distribution distribution_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.distribution ALTER COLUMN distribution_id SET DEFAULT nextval('public.distribution_distribution_id_seq'::regclass);


--
-- TOC entry 2706 (class 2604 OID 29803)
-- Name: protocol protocol_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.protocol ALTER COLUMN protocol_id SET DEFAULT nextval('public.protocol_protocol_id_seq'::regclass);


--
-- TOC entry 2707 (class 2604 OID 29804)
-- Name: survey survey_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.survey ALTER COLUMN survey_id SET DEFAULT nextval('public.survey_survey_id_seq'::regclass);


--
-- TOC entry 2714 (class 2606 OID 29806)
-- Name: dataset_keyword dataset_keyword_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.dataset_keyword
    ADD CONSTRAINT dataset_keyword_pkey PRIMARY KEY (dataset_id, keyword_id);


--
-- TOC entry 2709 (class 2606 OID 29858)
-- Name: dataset dataset_metadata_uuid_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.dataset
    ADD CONSTRAINT dataset_metadata_uuid_key UNIQUE (metadata_uuid);


--
-- TOC entry 2711 (class 2606 OID 29808)
-- Name: dataset dataset_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.dataset
    ADD CONSTRAINT dataset_pkey PRIMARY KEY (dataset_id);


--
-- TOC entry 2718 (class 2606 OID 29860)
-- Name: distribution distribution_distribution_url_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.distribution
    ADD CONSTRAINT distribution_distribution_url_key UNIQUE (distribution_url);


--
-- TOC entry 2720 (class 2606 OID 29810)
-- Name: distribution distribution_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.distribution
    ADD CONSTRAINT distribution_pkey PRIMARY KEY (dataset_id, distribution_id);


--
-- TOC entry 2724 (class 2606 OID 29812)
-- Name: keyword keyword_keyword_value_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.keyword
    ADD CONSTRAINT keyword_keyword_value_key UNIQUE (keyword_value);


--
-- TOC entry 2726 (class 2606 OID 29814)
-- Name: keyword keyword_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.keyword
    ADD CONSTRAINT keyword_pkey PRIMARY KEY (keyword_id);


--
-- TOC entry 2700 (class 2606 OID 29815)
-- Name: dataset latitude_max_constraint; Type: CHECK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE public.dataset
    ADD CONSTRAINT latitude_max_constraint CHECK (((('-90'::integer)::double precision <= latitude_max) AND (latitude_max <= (90)::double precision))) NOT VALID;


--
-- TOC entry 2900 (class 0 OID 0)
-- Dependencies: 2700
-- Name: CONSTRAINT latitude_max_constraint ON dataset; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON CONSTRAINT latitude_max_constraint ON public.dataset IS 'Range check constraint for maximum latitude';


--
-- TOC entry 2701 (class 2606 OID 29816)
-- Name: dataset latitude_min_constraint; Type: CHECK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE public.dataset
    ADD CONSTRAINT latitude_min_constraint CHECK (((('-90'::integer)::double precision <= latitude_min) AND (latitude_min <= (90)::double precision))) NOT VALID;


--
-- TOC entry 2901 (class 0 OID 0)
-- Dependencies: 2701
-- Name: CONSTRAINT latitude_min_constraint ON dataset; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON CONSTRAINT latitude_min_constraint ON public.dataset IS 'Range check constraint for minimum latitude';


--
-- TOC entry 2702 (class 2606 OID 29817)
-- Name: dataset longitude_max_constraint; Type: CHECK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE public.dataset
    ADD CONSTRAINT longitude_max_constraint CHECK (((('-180'::integer)::double precision < longitude_max) AND (longitude_max <= (180)::double precision))) NOT VALID;


--
-- TOC entry 2902 (class 0 OID 0)
-- Dependencies: 2702
-- Name: CONSTRAINT longitude_max_constraint ON dataset; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON CONSTRAINT longitude_max_constraint ON public.dataset IS 'Range check constraint for maximum longitude';


--
-- TOC entry 2703 (class 2606 OID 29818)
-- Name: dataset longitude_min_constraint; Type: CHECK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE public.dataset
    ADD CONSTRAINT longitude_min_constraint CHECK (((('-180'::integer)::double precision < longitude_min) AND (longitude_min <= (180)::double precision))) NOT VALID;


--
-- TOC entry 2903 (class 0 OID 0)
-- Dependencies: 2703
-- Name: CONSTRAINT longitude_min_constraint ON dataset; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON CONSTRAINT longitude_min_constraint ON public.dataset IS 'Range check constraint for minimum longitude';


--
-- TOC entry 2728 (class 2606 OID 29820)
-- Name: protocol protocol_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.protocol
    ADD CONSTRAINT protocol_pkey PRIMARY KEY (protocol_id);


--
-- TOC entry 2730 (class 2606 OID 29822)
-- Name: protocol protocol_protocol_value_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.protocol
    ADD CONSTRAINT protocol_protocol_value_key UNIQUE (protocol_value);


--
-- TOC entry 2732 (class 2606 OID 29824)
-- Name: survey survey_ga_survey_id_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.survey
    ADD CONSTRAINT survey_ga_survey_id_key UNIQUE (ga_survey_id);


--
-- TOC entry 2734 (class 2606 OID 29826)
-- Name: survey survey_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.survey
    ADD CONSTRAINT survey_pkey PRIMARY KEY (survey_id);


--
-- TOC entry 2715 (class 1259 OID 29827)
-- Name: fki_dataset_keyword_dataset_id; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX fki_dataset_keyword_dataset_id ON public.dataset_keyword USING btree (dataset_id);


--
-- TOC entry 2716 (class 1259 OID 29828)
-- Name: fki_dataset_keyword_keyword_id; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX fki_dataset_keyword_keyword_id ON public.dataset_keyword USING btree (keyword_id);


--
-- TOC entry 2712 (class 1259 OID 29829)
-- Name: fki_dataset_survey_id; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX fki_dataset_survey_id ON public.dataset USING btree (survey_id);


--
-- TOC entry 2721 (class 1259 OID 29830)
-- Name: fki_distribution_dataset_id; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX fki_distribution_dataset_id ON public.distribution USING btree (dataset_id);


--
-- TOC entry 2722 (class 1259 OID 29831)
-- Name: fki_distribution_protocol_id; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX fki_distribution_protocol_id ON public.distribution USING btree (protocol_id);


--
-- TOC entry 2736 (class 2606 OID 29832)
-- Name: dataset_keyword dataset_keyword_dataset_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.dataset_keyword
    ADD CONSTRAINT dataset_keyword_dataset_id_fkey FOREIGN KEY (dataset_id) REFERENCES public.dataset(dataset_id);


--
-- TOC entry 2737 (class 2606 OID 29837)
-- Name: dataset_keyword dataset_keyword_keyword_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.dataset_keyword
    ADD CONSTRAINT dataset_keyword_keyword_id_fkey FOREIGN KEY (keyword_id) REFERENCES public.keyword(keyword_id);


--
-- TOC entry 2735 (class 2606 OID 29842)
-- Name: dataset dataset_survey_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.dataset
    ADD CONSTRAINT dataset_survey_id_fkey FOREIGN KEY (survey_id) REFERENCES public.survey(survey_id);


--
-- TOC entry 2738 (class 2606 OID 29847)
-- Name: distribution distribution_dataset_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.distribution
    ADD CONSTRAINT distribution_dataset_id_fkey FOREIGN KEY (dataset_id) REFERENCES public.dataset(dataset_id);


--
-- TOC entry 2739 (class 2606 OID 29852)
-- Name: distribution distribution_protocol_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.distribution
    ADD CONSTRAINT distribution_protocol_id_fkey FOREIGN KEY (protocol_id) REFERENCES public.protocol(protocol_id) ON UPDATE CASCADE;


-- Completed on 2018-07-20 08:55:41

--
-- PostgreSQL database dump complete
--

